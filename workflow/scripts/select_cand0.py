#!/usr/bin/env python3
"""
Select the _0 TSS candidates (peaks below the spike-in threshold) to scan
for promoter motifs.

Inputs (from annotateTSS, tss-candidates branch)
  --peaks   {sample}_enriched_{ident}_peaks.bed   one window per candidate,
            name ends in _1 (passes threshold) or _0 (below threshold)
  --seqs    TSS_seq_{sample}_enriched_{ident}.fa.out   the same windows as
            sequences, headers  NAME::RECORD:start-end(strand), oriented 5'->3'
            along the transcript. meme_training.fa is the _1 subset of this
            file, so _0 windows taken from it have exactly the same geometry.
  --ratios  enr_ratios_{sample}_{ident}.plus.csv / .minus.csv (optional)
            per candidate: enriched read count (column 6) and enrichment
            ratio (last column). Only needed for the read-support filters.

Filters, in order (each is reported in the log)
  1. name ends in _0
  2. on the phage record(s) (--records)
  3. enriched reads >= --min-reads and enrichment ratio >= --min-ratio
  4. TSS not within --exclude-near-1 nt of a _1 TSS on the same strand
     (0 = keep all)

Outputs
  --out-bed    cand0.bed  (BED6, same columns as the peaks file)
  --out-fa     cand0.fa   (windows of the selected candidates, header = candidate name)
  --out-scan-fa  cand0_scan.fa  (same windows, header = short ID cand0001, ...)
  --out-ids    cand0_ids.tsv  (ID -> candidate name)
               FIMO 5.5 reads names such as "X::chr:1-55(+)" or "TSS_POS_4_0" as genomic
               coordinates and shifts the hit positions; plain IDs avoid that, and
               scan_candidates maps the IDs back to the candidate names afterwards.
  --out-table  cand0_table.tsv  (name, TSS, reads, ratio, filter outcome)
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict


def log(msg):
    print(f"[select_cand0] {msg}", file=sys.stderr)


def read_peaks(path):
    rows = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                sys.exit(f"[select_cand0] {path}:{n}: expected BED6, got {len(p)} columns")
            rows.append(dict(chrom=p[0], start=int(p[1]), end=int(p[2]), name=p[3],
                             strand=p[5].strip(), fields=p))
    return rows


def read_fasta(path):
    seqs, order, name, buf = {}, [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                # key by candidate name: "TSS_POS_2_1::LUZ7:57-112(+)" -> "TSS_POS_2_1"
                name, buf = line[1:].split()[0].split("::")[0], []
                order.append(name)
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs, order


def read_ratios(paths):
    """Rows: chrom, start, end, pos, strand, count, ... , ratio (last column)."""
    out = defaultdict(list)  # (chrom, strand) -> [(pos, reads, ratio)]
    for path in paths:
        with open(path, newline="") as fh:
            for row in csv.reader(fh):
                if len(row) < 7:
                    continue
                try:
                    chrom, pos, strand = row[0], int(float(row[3])), row[4].strip()
                    reads = float(row[5])
                    ratio = float(row[-1]) if row[-1] not in ("", "NA", "Inf", "inf") else float("inf")
                except ValueError:
                    continue  # header or malformed line
                out[(chrom, strand)].append((pos, reads, ratio))
    return out


def attach_ratios(peaks, ratios):
    """Find each candidate's ratio row: same record and strand, TSS inside the
    window, at the offset shared by most candidates on that strand."""
    offsets = {"+": Counter(), "-": Counter()}
    for p in peaks:
        for pos, _, _ in ratios.get((p["chrom"], p["strand"]), []):
            if p["start"] <= pos < p["end"]:
                offsets[p["strand"]][pos - p["start"]] += 1
    best = {s: (c.most_common(1)[0][0] if c else None) for s, c in offsets.items()}
    for s, off in best.items():
        if off is not None:
            log(f"strand {s}: TSS sits at window offset {off} "
                f"(window start + {off}; {offsets[s][off]} candidates agree)")
    index = {k: {pos: (r, q) for pos, r, q in v} for k, v in ratios.items()}
    for p in peaks:
        off = best.get(p["strand"])
        hit = None
        if off is not None:
            pos = p["start"] + off
            hit = index.get((p["chrom"], p["strand"]), {}).get(pos)
            p["tss"] = pos
        p["reads"], p["ratio"] = hit if hit else (None, None)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--seqs", required=True)
    ap.add_argument("--ratios", nargs="*", default=[])
    ap.add_argument("--records", nargs="+", required=True, help="phage record name(s)")
    ap.add_argument("--suffix", default="_0")
    ap.add_argument("--min-reads", type=float, default=0)
    ap.add_argument("--min-ratio", type=float, default=0)
    ap.add_argument("--exclude-near-1", type=int, default=0)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--out-fa", required=True)
    ap.add_argument("--out-scan-fa", required=True)
    ap.add_argument("--out-ids", required=True)
    ap.add_argument("--out-table", required=True)
    a = ap.parse_args()

    peaks = read_peaks(a.peaks)
    seqs, _ = read_fasta(a.seqs)
    records = set(a.records)
    present = {p["chrom"] for p in peaks}
    if not records & present:
        sys.exit(f"[select_cand0] none of the records {sorted(records)} occur in {a.peaks}; "
                 f"records present: {', '.join(sorted(present))}")

    have_ratios = bool(a.ratios)
    if have_ratios:
        attach_ratios(peaks, read_ratios(a.ratios))
    else:
        for p in peaks:
            p["reads"] = p["ratio"] = None
            p["tss"] = None
    if (a.min_reads or a.min_ratio) and not have_ratios:
        sys.exit("[select_cand0] --min-reads/--min-ratio need the enr_ratios csv files")

    ones = defaultdict(list)
    for p in peaks:
        if p["name"].endswith("_1") and p["chrom"] in records:
            ones[(p["chrom"], p["strand"])].append(p.get("tss") if p.get("tss") is not None
                                                   else (p["start"] + p["end"]) // 2)

    counts = Counter()
    kept = []
    with open(a.out_table, "w") as tab:
        tab.write("name\tchrom\twindow_start\twindow_end\tstrand\ttss\tenriched_reads\t"
                  "enrichment_ratio\toutcome\n")
        for p in peaks:
            if not p["name"].endswith(a.suffix):
                continue
            counts["_0 candidates (all records)"] += 1
            if p["chrom"] not in records:
                continue
            counts["on phage record(s)"] += 1
            if p["name"] not in seqs:
                outcome = "no_sequence_in_fa_out"
            elif have_ratios and p["reads"] is None and (a.min_reads or a.min_ratio):
                outcome = "no_ratio_row"
            elif a.min_reads and p["reads"] < a.min_reads:
                outcome = "low_reads"
            elif a.min_ratio and p["ratio"] < a.min_ratio:
                outcome = "low_ratio"
            else:
                tss = p.get("tss") if p.get("tss") is not None else (p["start"] + p["end"]) // 2
                near = [t for t in ones[(p["chrom"], p["strand"])] if abs(t - tss) <= a.exclude_near_1]
                outcome = "near_1" if (a.exclude_near_1 and near) else "selected"
            counts[outcome] += 1
            if outcome == "selected":
                kept.append(p)
            fmt = lambda v: "" if v is None else (f"{v:g}" if isinstance(v, float) else str(v))
            tab.write(f"{p['name']}\t{p['chrom']}\t{p['start']}\t{p['end']}\t{p['strand']}\t"
                      f"{fmt(p.get('tss'))}\t{fmt(p['reads'])}\t{fmt(p['ratio'])}\t{outcome}\n")

    with open(a.out_bed, "w") as bed, open(a.out_fa, "w") as fa, \
            open(a.out_scan_fa, "w") as sfa, open(a.out_ids, "w") as ids:
        ids.write("id\tname\n")
        for i, p in enumerate(kept, 1):
            bed.write("\t".join(p["fields"]) + "\n")
            fa.write(f">{p['name']}\n{seqs[p['name']]}\n")
            # short ID without separators or number pairs, so FIMO cannot read it
            # as a genomic coordinate (see --out-scan-fa above)
            sid = f"cand{i:04d}"
            sfa.write(f">{sid}\n{seqs[p['name']]}\n")
            ids.write(f"{sid}\t{p['name']}\n")

    for k in ["_0 candidates (all records)", "on phage record(s)", "no_sequence_in_fa_out",
              "no_ratio_row", "low_reads", "low_ratio", "near_1", "selected"]:
        if counts[k] or k in ("on phage record(s)", "selected"):
            log(f"{k:28s} {counts[k]}")
    if counts["no_sequence_in_fa_out"]:
        log("WARNING: some _0 candidates have no window in the .fa.out file")
    if not kept:
        log("WARNING: no _0 candidates selected; cand0.fa is empty")


if __name__ == "__main__":
    main()
