#!/usr/bin/env python3
"""
rescue_tss: rescue _0 TSS candidates that have a promoter motif at the right place.

0. Motif selection (motif_selection.tsv)
   For every motif in the motif file, MEME's training sites (the _1 windows) give
   the motif's E-value and where it sits in the window. A motif is used when
     - it has >= --min-training-sites sites,
     - its sites are positionally conserved: max start - min start <= --max-site-spread
       (a real promoter element sits at a fixed distance from the TSS),
     - its E-value <= --max-evalue (if given).
   --motifs overrides the automatic choice (e.g. --motifs MEME-1).
   Position rule: a hit only counts if its start lies in
     [min training start - slack, max training start + slack]   (--position auto)
   or in a fixed range (--position 26 41), or anywhere (--position off).
   All positions are 1-based positions inside the window. _1, _0 and random
   windows have the same length and orientation, so one rule fits all three.

1. Calibration on the random windows (background_rates.tsv, calibrated_cutoffs.tsv)
   Hits are filtered with the SAME position rule as the candidates, so the
   background rate is that of the rule actually applied.
   a. per motif: lower the p-value cut-off (start = FIMO --thresh) by /1.1 until at
      most --target-fpr of the random windows have a hit;
   b. jointly: while the fraction of random windows with a hit of ANY used motif is
      above --target-fpr, lower all cut-offs by /1.1.
   So --target-fpr is the false-positive rate of the rescue as a whole.

2. Rescue: a _0 candidate is rescued if any used motif hits it inside its position
   range with p <= that motif's cut-off.

3. Tiers: 1 = _1 TSS (above the spike-in threshold), 2 = rescued _0, 3 = _0 not rescued.
"""
import argparse
import os
import re
import sys

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--random-fimo", required=True)   # FIMO hits in random windows
parser.add_argument("--random-fa", required=True)     # all random windows
parser.add_argument("--cand0-fimo", required=True)    # FIMO hits in _0 candidates
parser.add_argument("--cand0-bed", required=True)     # scanned _0 candidates
parser.add_argument("--tss-fa", required=True)        # all TSS windows (_1 and _0)
parser.add_argument("--motif-file", required=True)    # meme.txt (training sites, E-values)
parser.add_argument("--records", nargs="*", default=None)  # phage records; default: those in cand0.bed
parser.add_argument("--target-fpr", type=float, default=0.05)
parser.add_argument("--start-cutoff", type=float, default=0.001)  # = FIMO --thresh
parser.add_argument("--motifs", nargs="*", default=None)          # override selection
parser.add_argument("--max-evalue", type=float, default=None)
parser.add_argument("--min-training-sites", type=int, default=3)
parser.add_argument("--max-site-spread", type=int, default=12)
parser.add_argument("--position", nargs="+", default=["auto"])   # auto | off | MIN MAX
parser.add_argument("--position-slack", type=int, default=3)
parser.add_argument("--up", type=int, default=None)    # TSS sequence extraction: upstream
parser.add_argument("--down", type=int, default=None)  # TSS sequence extraction: downstream
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)


def read_fimo(path):
    d = pd.read_csv(path, sep="\t")
    d = d[~d["motif_id"].astype(str).str.startswith("#")]   # FIMO comment lines, if any
    d["start"] = d["start"].astype(int)
    d["p-value"] = d["p-value"].astype(float)
    return d


# ---------------------------------------------------------------------------
# 0. Motifs: E-value and training-site positions from meme.txt
# ---------------------------------------------------------------------------
def read_meme_sites(path):
    """{alt_id: {"id", "evalue", "sites": [1-based window starts]}} from meme.txt.
    Motif files without a MEME 'sites sorted by position' section give no sites."""
    motifs, current, in_sites = {}, None, False
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^MOTIF\s+(\S+)\s+(\S+)?(.*)$", line)
            if m:
                mid, alt = m.group(1), m.group(2) or m.group(1)
                ev = re.search(r"E-value\s*=\s*(\S+)", m.group(3) or "")
                motifs[alt] = {"id": mid, "evalue": float(ev.group(1)) if ev else None, "sites": []}
                current, in_sites = alt, False
                continue
            if "sites sorted by position p-value" in line:
                in_sites = None          # header follows: "Sequence name ...", then a dash line
                continue
            if in_sites is None and line.startswith("Sequence name"):
                in_sites = "header"
                continue
            if in_sites == "header":     # the "-------------  -----" line under the header
                in_sites = True
                continue
            if in_sites:
                if line.startswith("----") or not line.strip():
                    in_sites = False
                    continue
                tok = line.split()
                # name [strand] start p-value ...  (strand only with -revcomp)
                start = tok[2] if tok[1] in ("+", "-") else tok[1]
                motifs[current]["sites"].append(int(start))
    return motifs


meme = read_meme_sites(args.motif_file)
rnd = read_fimo(args.random_fimo)
cand = read_fimo(args.cand0_fimo)

# motifs FIMO knows about (a motif file without a MOTIF block gives none)
all_motifs = list(dict.fromkeys(list(meme) + list(rnd["motif_alt_id"].unique())))

sel_rows, ranges = [], {}
for alt in all_motifs:
    info = meme.get(alt, {"id": alt, "evalue": None, "sites": []})
    sites = info["sites"]
    spread = (max(sites) - min(sites)) if sites else None

    # position range
    if args.position[0] == "off":
        rng = None
    elif args.position[0] == "auto":
        rng = (min(sites) - args.position_slack, max(sites) + args.position_slack) if sites else None
    else:
        rng = (int(args.position[0]), int(args.position[1]))

    # selection
    if args.motifs:
        use, why = alt in args.motifs, "listed in config" if alt in args.motifs else "not listed in config"
    elif not sites:
        use, why = True, "no MEME training sites in motif file: used without selection"
    elif len(sites) < args.min_training_sites:
        use, why = False, f"only {len(sites)} training sites (< {args.min_training_sites})"
    elif spread > args.max_site_spread:
        use, why = False, f"not positional: training sites spread over {spread} nt (> {args.max_site_spread})"
    elif args.max_evalue is not None and info["evalue"] is not None and info["evalue"] > args.max_evalue:
        use, why = False, f"E-value {info['evalue']:.2g} > {args.max_evalue}"
    else:
        use, why = True, "positional"
    if use:
        ranges[alt] = rng

    sel_rows.append(dict(
        motif=alt, consensus=info["id"], evalue=info["evalue"], n_training_sites=len(sites),
        training_starts=",".join(map(str, sorted(sites))), spread=spread,
        allowed_start_min=rng[0] if rng else None, allowed_start_max=rng[1] if rng else None,
        used=use, reason=why))
    print(f"{alt}: {'USED' if use else 'dropped'} ({why}); allowed starts {rng}")

sel = pd.DataFrame(sel_rows)
sel.to_csv(f"{args.outdir}/motif_selection.tsv", sep="\t", index=False)
if not ranges:
    print("WARNING: no motif selected -> nothing can be rescued")
if all(e is not None and e > 0.05 for e in sel.loc[sel["used"], "evalue"]) and ranges:
    print("NOTE: no used motif is significant in MEME (E > 0.05); pool more _1 windows")


def apply_rules(d):
    """keep hits of used motifs inside their position range"""
    d = d[d["motif_alt_id"].isin(ranges)].copy()
    lo = d["motif_alt_id"].map(lambda m: ranges[m][0] if ranges[m] else -10**9)
    hi = d["motif_alt_id"].map(lambda m: ranges[m][1] if ranges[m] else 10**9)
    return d[(d["start"] >= lo) & (d["start"] <= hi)]


rnd_ok = apply_rules(rnd)
cand_ok = apply_rules(cand)


# ---------------------------------------------------------------------------
# 1. Calibration on the random windows
# ---------------------------------------------------------------------------
n_windows = 0
with open(args.random_fa) as f:
    for line in f:
        if line.startswith(">"):
            n_windows += 1
print(f"random windows: {n_windows}")
if n_windows == 0:
    sys.exit("no random windows in " + args.random_fa)


def rate(hits):
    return hits["sequence_name"].nunique() / n_windows


def combined_rate(cuts):
    return rate(rnd_ok[rnd_ok["p-value"] <= rnd_ok["motif_alt_id"].map(cuts)])


rows = []
cutoffs = {}
for motif in ranges:                                   # a. per motif
    m = rnd_ok[rnd_ok["motif_alt_id"] == motif]
    cutoff = args.start_cutoff
    while True:
        fpr = rate(m[m["p-value"] <= cutoff])
        rows.append((motif, "per motif", cutoff, int(fpr * n_windows), fpr))
        if fpr <= args.target_fpr:
            break
        cutoff /= 1.1
    cutoffs[motif] = cutoff

combined = combined_rate(cutoffs) if cutoffs else 0.0     # b. jointly
while cutoffs and combined > args.target_fpr:
    cutoffs = {k: v / 1.1 for k, v in cutoffs.items()}
    combined = combined_rate(cutoffs)
    rows.append(("any", "joint", None, int(combined * n_windows), combined))
for motif, c in cutoffs.items():
    print(f"{motif}: cutoff {c:.3g}, background rate {rate(rnd_ok[(rnd_ok['motif_alt_id'] == motif) & (rnd_ok['p-value'] <= c)]):.3f}")
print(f"combined background rate (any motif): {combined:.3f}")


# ---------------------------------------------------------------------------
# 2. Rescue the _0 candidates
# ---------------------------------------------------------------------------
cand_ok = cand_ok.copy()
cand_ok["cutoff"] = cand_ok["motif_alt_id"].map(cutoffs)
passed = cand_ok[cand_ok["p-value"] <= cand_ok["cutoff"]]
best = passed.sort_values("p-value").drop_duplicates("sequence_name").set_index("sequence_name")

bed0 = pd.read_csv(args.cand0_bed, sep="\t", header=None, comment="#")
bed0 = bed0.iloc[:, :6]
bed0.columns = ["chrom", "start", "end", "name", "score", "strand"]

if len(cand) > 0 and not set(cand["sequence_name"]) & set(bed0["name"]):
    print("WARNING: no _0 candidate name in the BED file matches a FIMO "
          "sequence name; check the naming of cand0.bed and fimo_cand0.tsv")

bed0["tier"] = bed0["name"].apply(lambda n: 2 if n in best.index else 3)
bed0["best_motif"] = bed0["name"].map(best["motif_alt_id"])
bed0["best_p"] = bed0["name"].map(best["p-value"])
bed0["motif_start_in_window"] = bed0["name"].map(best["start"])


# ---------------------------------------------------------------------------
# 3. The _1 TSSs from the TSS FASTA (tier 1)
# ---------------------------------------------------------------------------
# header example: >TSS_POS_38_1::LUZ7:8210-8265(+)
keep_records = set(args.records) if args.records else set(bed0["chrom"].astype(str))
print(f"records kept for tier 1: {sorted(keep_records)}")

tss1 = []
with open(args.tss_fa) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        name, location = line[1:].strip().split("::")
        if name.split("_")[-1] != "1":
            continue
        chrom, coords = location.rsplit(":", 1)
        if chrom not in keep_records:
            continue
        strand = coords[-2]
        start, end = coords[:-3].split("-")
        tss1.append((chrom, int(start), int(end), name, 0, strand))

bed1 = pd.DataFrame(tss1, columns=["chrom", "start", "end", "name", "score", "strand"])
bed1["tier"] = 1
bed1["best_motif"] = pd.NA
bed1["best_p"] = pd.NA
bed1["motif_start_in_window"] = pd.NA


# ---------------------------------------------------------------------------
# 4. Output files
# ---------------------------------------------------------------------------
final = pd.concat([bed1, bed0], ignore_index=True)
# TSS position, same convention as the enr_ratios csv (window = TSS-up .. TSS+down)
if args.up is not None and args.down is not None:
    final["tss"] = final.apply(
        lambda r: r["start"] + (args.up if r["strand"] == "+" else args.down), axis=1)

final.to_csv(f"{args.outdir}/tss_tiers.tsv", sep="\t", index=False)
final[final["tier"] <= 2][["chrom", "start", "end", "name", "score", "strand"]].to_csv(
    f"{args.outdir}/tss_final.bed", sep="\t", index=False, header=False)
pd.DataFrame(rows, columns=["motif", "step", "p_cutoff", "n_hit_windows", "fpr"]).to_csv(
    f"{args.outdir}/background_rates.tsv", sep="\t", index=False)
pd.DataFrame(list(cutoffs.items()), columns=["motif", "p_cutoff"]).to_csv(
    f"{args.outdir}/calibrated_cutoffs.tsv", sep="\t", index=False)

n0 = len(bed0)
n_rescued = int((bed0["tier"] == 2).sum())
expected_false = combined * n0
with open(f"{args.outdir}/rescue_summary.txt", "w") as f:
    f.write(f"motifs used\t{','.join(ranges) or 'none'}\n")
    for m, r in ranges.items():
        f.write(f"allowed starts {m}\t{r[0]}-{r[1]}\n" if r else f"allowed starts {m}\tany\n")
    f.write(f"random windows\t{n_windows}\n")
    f.write(f"target fpr (any motif)\t{args.target_fpr}\n")
    f.write(f"combined background rate\t{combined:.4f}\n")
    f.write(f"_1 TSSs (tier 1)\t{len(bed1)}\n")
    f.write(f"_0 candidates scanned\t{n0}\n")
    f.write(f"_0 rescued (tier 2)\t{n_rescued}\n")
    f.write(f"expected false rescues\t{expected_false:.1f}\n")
    if n_rescued > 0:
        f.write(f"estimated FDR among rescued\t{expected_false / n_rescued:.3f}\n")

print(open(f"{args.outdir}/rescue_summary.txt").read())
