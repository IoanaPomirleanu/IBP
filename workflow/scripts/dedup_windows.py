#!/usr/bin/env python3
"""
Collapse redundant promoter windows before motif discovery.

Input headers are as written by `bedtools getfasta -name -s`:
    >TSS_POS_2_1::LUZ7:57-112(+)

Two windows are redundant if
  - they are on the same record and strand and their starts are within
    --max-offset nt of each other (the same promoter called twice), or
  - their sequences are identical (e.g. copies of a promoter in the
    terminal repeats at both genome ends).
The first window of each redundant group (lowest coordinate) is kept.

MEME treats every sequence as independent evidence, so duplicates make a
repeated sequence look like a conserved motif.
"""
import argparse
import re
import sys

HEADER = re.compile(r"^(?P<name>.+)::(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)$")


def read_fasta(path):
    recs, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    recs.append((name, "".join(seq).upper()))
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name:
        recs.append((name, "".join(seq).upper()))
    return recs


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--max-offset", type=int, default=15,
                   help="windows starting within this many nt are one promoter (use the TSS cluster width)")
    a = p.parse_args()

    recs = []
    for header, seq in read_fasta(a.fasta):
        m = HEADER.match(header)
        if not m:
            sys.exit(f"[dedup_windows] unexpected header '{header}' "
                     f"(expected name::chrom:start-end(strand))")
        recs.append(dict(id=header, seq=seq, chrom=m["chrom"], strand=m["strand"],
                         start=int(m["start"])))
    recs.sort(key=lambda r: (r["chrom"], r["strand"], r["start"]))

    # 1. group windows of the same promoter (same record/strand, starts close together)
    groups = []
    for r in recs:
        g = groups[-1] if groups else None
        if g and g[0]["chrom"] == r["chrom"] and g[0]["strand"] == r["strand"] \
                and r["start"] - g[-1]["start"] <= a.max_offset:
            g.append(r)
        else:
            groups.append([r])

    # 2. drop groups that share any sequence with an already kept group (repeat copies)
    kept, kept_seqs, log = [], {}, []
    for g in groups:
        rep = g[0]
        for r in g[1:]:
            log.append(f"removed {r['id']}: within {a.max_offset} nt of {rep['id']}")
        shared = next((kept_seqs[r["seq"]] for r in g if r["seq"] in kept_seqs), None)
        if shared:
            log.append(f"removed {rep['id']}: same promoter as {shared} (repeat copy)")
            continue
        kept.append(rep)
        for r in g:
            kept_seqs[r["seq"]] = rep["id"]

    with open(a.out, "w") as out:
        for r in kept:
            out.write(f">{r['id']}\n{r['seq']}\n")
    for line in log:
        print(f"[dedup_windows] {line}", file=sys.stderr)
    print(f"[dedup_windows] kept {len(kept)} of {len(recs)} windows", file=sys.stderr)


if __name__ == "__main__":
    main()
