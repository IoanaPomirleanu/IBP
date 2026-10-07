#!/usr/bin/env python3
"""
Write only selected records of a FASTA file to a new FASTA.

Used to build the motif background model from the phage alone when the
pipeline's genome FASTA is a concatenation (phage + host + spike-in).

usage: extract_records.py IN.fa OUT.fa [RECORD_ID ...]
With no record IDs, every record is copied.
"""
import sys


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    src, out, keep = sys.argv[1], sys.argv[2], set(sys.argv[3:])
    found, all_ids, write = set(), [], False
    with open(src) as fi, open(out, "w") as fo:
        for line in fi:
            if line.startswith(">"):
                rid = line[1:].split()[0]
                all_ids.append(rid)
                write = not keep or rid in keep
                if write:
                    found.add(rid)
            if write:
                fo.write(line)
    missing = keep - found
    if missing:
        sys.exit(f"[extract_records] record(s) {sorted(missing)} not in {src}; "
                 f"records present: {', '.join(all_ids)}")
    print(f"[extract_records] background genome records: {', '.join(sorted(found))}")


if __name__ == "__main__":
    main()
