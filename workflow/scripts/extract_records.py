#!/usr/bin/env python3
"""
Write only selected records of a FASTA file to a new FASTA.

Used to build the motif background model from the phage alone when the
pipeline's genome FASTA is a concatenation (phage + host + spike-in).

  extract_records.py IN.fa OUT.fa --records LUZ7 [...]
  extract_records.py IN.fa OUT.fa --records-from meme_training.fa
      (record names taken from headers NAME::RECORD:start-end(strand))
With neither option, every record is copied.
"""
import argparse
import sys


def records_from_headers(path):
    recs = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">") and "::" in line:
                recs.add(line[1:].split()[0].split("::", 1)[1].rsplit(":", 1)[0])
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("out")
    ap.add_argument("--records", nargs="*", default=[])
    ap.add_argument("--records-from", default=None)
    a = ap.parse_args()

    keep = set(a.records)
    if not keep and a.records_from:
        keep = records_from_headers(a.records_from)
        if not keep:
            sys.exit(f"[extract_records] no 'NAME::RECORD:start-end(strand)' headers in {a.records_from}")

    found, all_ids, write = set(), [], False
    with open(a.src) as fi, open(a.out, "w") as fo:
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
        sys.exit(f"[extract_records] record(s) {sorted(missing)} not in {a.src}; "
                 f"records present: {', '.join(all_ids)}. "
                 "Set 'genome fasta' / 'genome records' in the motifs: config block.")
    print(f"[extract_records] background genome records: {', '.join(sorted(found))}")


if __name__ == "__main__":
    main()
