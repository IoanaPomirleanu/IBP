#!/usr/bin/env python3
"""
rescue_tss: rescue _0 TSS candidates that have a promoter motif.

1. Calibration: for each motif, lower the p-value cutoff (starting at the
   FIMO threshold 0.001) until at most --target-fpr of the random genomic
   windows still have a hit.
2. Rescue: a _0 candidate is rescued if any motif hits it with a p-value
   at or below that motif's cutoff.
3. Tiers: 1 = _1 TSS (above the original threshold), 2 = rescued _0,
   3 = _0 not rescued.
"""
import argparse
import os

import pandas as pd

# ---------------------------------------------------------------------------
# Options from the Snakemake rule
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--random-fimo", required=True)   # FIMO hits in random windows
parser.add_argument("--random-fa", required=True)     # all random windows
parser.add_argument("--cand0-fimo", required=True)    # FIMO hits in _0 candidates
parser.add_argument("--cand0-bed", required=True)     # scanned _0 candidates
parser.add_argument("--tss-fa", required=True)        # all TSS windows (_1 and _0)
parser.add_argument("--records", nargs="*", default=None)  # phage records; default: those in cand0.bed
parser.add_argument("--target-fpr", type=float, default=0.05)
parser.add_argument("--start-cutoff", type=float, default=0.001)  # = FIMO --thresh
parser.add_argument("--outdir", required=True)
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Calibration on the random windows
# ---------------------------------------------------------------------------
# number of random windows = number of headers in random.fa
n_windows = 0
with open(args.random_fa) as f:
    for line in f:
        if line.startswith(">"):
            n_windows += 1
print(f"random windows: {n_windows}")

df = pd.read_csv(args.random_fimo, sep="\t")
df = df[~df["motif_id"].astype(str).str.startswith("#")]   # skip FIMO comment lines, if any

cutoffs = {}   # final cutoff per motif
rows = []      # every step, for background_rates.tsv

for motif in df["motif_alt_id"].unique():
    m = df[df["motif_alt_id"] == motif]       # only this motif's hits
    cutoff = args.start_cutoff                # start again for each motif

    while True:
        hits = m[m["p-value"] <= cutoff]
        n_hit_windows = hits["sequence_name"].nunique()   # windows, not hits
        fpr = n_hit_windows / n_windows
        rows.append((motif, cutoff, n_hit_windows, fpr))

        if fpr <= args.target_fpr:
            break
        cutoff = cutoff / 1.1

    cutoffs[motif] = cutoff
    print(f"{motif}: cutoff {cutoff:.3g}, background rate {fpr:.3f}")

# Combined background rate: fraction of random windows where ANY motif passes
# its own cutoff. This is the real false-positive rate of the rescue, because
# a candidate is rescued as soon as one motif passes.
bg_pass = df[df["p-value"] <= df["motif_alt_id"].map(cutoffs)]
combined_fpr = bg_pass["sequence_name"].nunique() / n_windows
print(f"combined background rate (any motif): {combined_fpr:.3f}")


# ---------------------------------------------------------------------------
# 2. Rescue the _0 candidates
# ---------------------------------------------------------------------------
cand = pd.read_csv(args.cand0_fimo, sep="\t")
cand = cand[~cand["motif_id"].astype(str).str.startswith("#")]

# keep only hits that pass their motif's cutoff (motifs without a cutoff never pass)
cand["cutoff"] = cand["motif_alt_id"].map(cutoffs)
passed = cand[cand["p-value"] <= cand["cutoff"]]

# best (lowest p-value) passing hit per candidate
best = passed.sort_values("p-value").drop_duplicates("sequence_name")
best = best.set_index("sequence_name")

# the scanned _0 candidates (BED: chrom, start, end, name, score, strand)
bed0 = pd.read_csv(args.cand0_bed, sep="\t", header=None, comment="#")
bed0 = bed0.iloc[:, :6]
bed0.columns = ["chrom", "start", "end", "name", "score", "strand"]

# check that names in the FIMO output and the BED file match
if not set(cand["sequence_name"]) & set(bed0["name"]) and len(cand) > 0:
    print("WARNING: no _0 candidate name in the BED file matches a FIMO "
          "sequence name; check the naming of cand0.bed and fimo_cand0.tsv")

bed0["tier"] = bed0["name"].apply(lambda n: 2 if n in best.index else 3)
bed0["best_motif"] = bed0["name"].map(best["motif_alt_id"])
bed0["best_p"] = bed0["name"].map(best["p-value"])


# ---------------------------------------------------------------------------
# 3. The _1 TSSs from the TSS FASTA (tier 1)
# ---------------------------------------------------------------------------
# header example: >TSS_POS_38_1::LUZ7:8210-8265(+)
# Keep only _1 TSSs on the phage: the same records as the scanned _0
# candidates (cand0.bed only contains phage records), unless --records is given.
keep_records = set(args.records) if args.records else set(bed0["chrom"].astype(str))
print(f"records kept for tier 1: {sorted(keep_records)}")

tss1 = []
with open(args.tss_fa) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        header = line[1:].strip()
        name, location = header.split("::")
        if name.split("_")[-1] != "1":          # keep only _1
            continue
        chrom, coords = location.rsplit(":", 1)
        if chrom not in keep_records:                # skip host TSSs
            continue
        strand = coords[-2]                     # "+" or "-"
        start, end = coords[:-3].split("-")
        tss1.append((chrom, int(start), int(end), name, 0, strand))

bed1 = pd.DataFrame(tss1, columns=["chrom", "start", "end", "name", "score", "strand"])
bed1["tier"] = 1
bed1["best_motif"] = pd.NA
bed1["best_p"] = pd.NA


# ---------------------------------------------------------------------------
# 4. Output files
# ---------------------------------------------------------------------------
final = pd.concat([bed1, bed0], ignore_index=True)

# all TSSs with their tier
final.to_csv(f"{args.outdir}/tss_tiers.tsv", sep="\t", index=False)

# final TSS set: tiers 1 and 2, as BED
final[final["tier"] <= 2][["chrom", "start", "end", "name", "score", "strand"]].to_csv(
    f"{args.outdir}/tss_final.bed", sep="\t", index=False, header=False)

# every calibration step
pd.DataFrame(rows, columns=["motif", "p_cutoff", "n_hit_windows", "fpr"]).to_csv(
    f"{args.outdir}/background_rates.tsv", sep="\t", index=False)

# final cutoff per motif
pd.DataFrame(list(cutoffs.items()), columns=["motif", "p_cutoff"]).to_csv(
    f"{args.outdir}/calibrated_cutoffs.tsv", sep="\t", index=False)

# summary: the numbers to report
n0 = len(bed0)
n_rescued = int((bed0["tier"] == 2).sum())
expected_false = combined_fpr * n0
with open(f"{args.outdir}/rescue_summary.txt", "w") as f:
    f.write(f"random windows\t{n_windows}\n")
    f.write(f"target fpr per motif\t{args.target_fpr}\n")
    f.write(f"combined background rate\t{combined_fpr:.4f}\n")
    f.write(f"_1 TSSs (tier 1)\t{len(bed1)}\n")
    f.write(f"_0 candidates scanned\t{n0}\n")
    f.write(f"_0 rescued (tier 2)\t{n_rescued}\n")
    f.write(f"expected false rescues\t{expected_false:.1f}\n")
    if n_rescued > 0:
        f.write(f"estimated FDR among rescued\t{expected_false / n_rescued:.3f}\n")

print(open(f"{args.outdir}/rescue_summary.txt").read())