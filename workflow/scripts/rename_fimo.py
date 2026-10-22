#!/usr/bin/env python3
"""
Replace the short window IDs in FIMO --text output (cand0001, ...) with the
candidate names from cand0_ids.tsv. Positions are left untouched: they are
1-based positions inside the window.

usage: rename_fimo.py FIMO_IN IDS_TSV FIMO_OUT
"""
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    fimo_in, ids_tsv, fimo_out = sys.argv[1:]
    names = {}
    with open(ids_tsv) as fh:
        next(fh, None)
        for line in fh:
            sid, name = line.rstrip("\n").split("\t")[:2]
            names[sid] = name

    col, n, unknown = None, 0, set()
    with open(fimo_in) as fi, open(fimo_out, "w") as fo:
        for line in fi:
            if col is None and line.strip():
                head = [h.strip().lstrip("#").strip() for h in line.rstrip("\n").split("\t")]
                if "sequence_name" not in head:
                    sys.exit(f"[rename_fimo] no 'sequence_name' column in {fimo_in}")
                col = head.index("sequence_name")
                fo.write(line)
                continue
            if not line.strip() or line.startswith("#"):
                fo.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if f[col] in names:
                f[col] = names[f[col]]
            else:
                unknown.add(f[col])
            fo.write("\t".join(f) + "\n")
            n += 1
    if unknown:
        sys.exit(f"[rename_fimo] FIMO reported unknown window IDs {sorted(unknown)[:5]}: "
                 "FIMO may have altered the sequence names")
    print(f"[rename_fimo] {n} hits, window IDs replaced by candidate names", file=sys.stderr)


if __name__ == "__main__":
    main()
