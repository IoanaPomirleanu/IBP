#!/usr/bin/env python3
"""
prepare_known_motifs: one MEME-format motif file for both FIMO scans, with
  - the de novo motifs from MEME (meme.txt, optional), unchanged (MEME-1, ...)
  - known motifs supplied by the user, renamed KNOWN-<name>:
      * from a MEME-format file (minimal MEME format, or MEME/STREME output;
        JASPAR / TRANSFAC files can be converted first with the MEME suite tools
        jaspar2meme / transfac2meme)
      * from IUPAC consensus sequences, e.g. sigma70_-10=TATAAT
        (each allowed base shares the probability; every other base gets
        --mismatch so that a site with a mismatch can still score)

Outputs
  --out        scan_motifs.meme   minimal MEME format (what FIMO reads)
  --out-list   known_motifs.tsv   alt id, name, width, source of every known motif
"""
import argparse
import re
import sys

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
         "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
BASES = "ACGT"
PREFIX = "KNOWN-"


def read_meme_motifs(path):
    """[(id, alt, rows)] from any MEME-format file (minimal or full MEME/STREME
    output): MOTIF line, then the first 'letter-probability matrix' after it."""
    motifs, cur, want, rows, w = [], None, False, [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("MOTIF"):
                p = line.split()
                cur = [p[1], p[2] if len(p) > 2 and not p[2].startswith("width") else "", None]
                motifs.append(cur)
                want = True
                continue
            if cur is not None and want and line.startswith("letter-probability matrix"):
                m = re.search(r"w=\s*(\d+)", line)
                w, rows, want = int(m.group(1)) if m else 0, [], False
                cur[2] = rows
                continue
            if cur is not None and cur[2] is rows and w and len(rows) < w:
                v = line.split()
                if len(v) == 4:
                    rows.append([float(x) for x in v])
    out = [(i, a, r) for i, a, r in motifs if r]
    for i, a, r in motifs:
        if not r:
            print(f"[prepare_known_motifs] {path}: motif {i} has no letter-probability matrix, skipped",
                  file=sys.stderr)
    return out


def consensus_rows(seq, mismatch):
    rows = []
    for ch in seq.upper():
        if ch not in IUPAC:
            sys.exit(f"[prepare_known_motifs] '{seq}': '{ch}' is not an IUPAC DNA letter")
        ok = IUPAC[ch]
        if len(ok) == 4:
            rows.append([0.25] * 4)
            continue
        miss = mismatch
        hit = (1 - miss * (4 - len(ok))) / len(ok)
        rows.append([hit if b in ok else miss for b in BASES])
    return rows


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.+-]", "_", name)


ap = argparse.ArgumentParser()
ap.add_argument("--denovo", default=None, help="meme.txt of discoverMotifs (optional)")
ap.add_argument("--known-file", default=None, help="MEME-format file with known motifs")
ap.add_argument("--consensus", nargs="*", default=[], help="NAME=IUPAC ...")
ap.add_argument("--mismatch", type=float, default=0.05)
ap.add_argument("--out", required=True)
ap.add_argument("--out-list", required=True)
a = ap.parse_args()

blocks = []                      # (id, alt, rows)
if a.denovo:
    dn = read_meme_motifs(a.denovo)
    print(f"[prepare_known_motifs] de novo motifs: {len(dn)}", file=sys.stderr)
    blocks += [(i, alt or i, r) for i, alt, r in dn]

known = []                       # (alt, name, width, source)
if a.known_file:
    for i, alt, r in read_meme_motifs(a.known_file):
        name = safe(alt or i)
        blocks.append((safe(i), PREFIX + name, r))
        known.append((PREFIX + name, alt or i, len(r), a.known_file))
for item in a.consensus:
    if "=" not in item:
        sys.exit(f"[prepare_known_motifs] consensus '{item}': expected NAME=SEQUENCE")
    name, seq = item.split("=", 1)
    r = consensus_rows(seq, a.mismatch)
    blocks.append((seq.upper(), PREFIX + safe(name), r))
    known.append((PREFIX + safe(name), name, len(r), "consensus " + seq.upper()))

alts = [b[1] for b in blocks]
dups = sorted({x for x in alts if alts.count(x) > 1})
if dups:
    sys.exit(f"[prepare_known_motifs] duplicate motif names: {', '.join(dups)}")
print(f"[prepare_known_motifs] known motifs: {len(known)}", file=sys.stderr)

with open(a.out, "w") as fh:
    fh.write("MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
             "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
    for i, alt, r in blocks:
        fh.write(f"MOTIF {i} {alt}\n\n")
        fh.write(f"letter-probability matrix: alength= 4 w= {len(r)} nsites= 20 E= 0\n")
        for row in r:
            s = sum(row)
            fh.write(" " + "  ".join(f"{x / s:.6f}" for x in row) + "\n")
        fh.write("\n")

with open(a.out_list, "w") as fh:
    fh.write("motif\tname\twidth\tsource\n")
    for alt, name, w, src in known:
        fh.write(f"{alt}\t{name}\t{w}\t{src}\n")
