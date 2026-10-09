# Promoter motif module: discovery and _0 candidate scan

## Rules

| File | Rule | What it does |
|---|---|---|
| `rules/motifs.smk` | `dedup_training` | collapses redundant `_1` windows (same promoter called twice, terminal-repeat copies) |
| | `motif_background_model` | order-2 background of the phage genome only (for MEME) |
| | `discover_motifs` | MEME (zoops, given strand only) on the deduplicated `_1` windows → `meme/meme.txt` |
| `rules/scan_candidates.smk` | `select_cand0` | picks the `_0` candidates on the phage record from the peaks BED and writes their windows (taken from the same `TSS_seq_*.fa.out` as `meme_training.fa`) |
| | `scan_candidates` | FIMO on the `_0` windows with the shared motif file, `bg.txt`, threshold and options from `motif background` (identical to `scan_background`), then window IDs → candidate names |

Everything runs from `rule all` (target `results/promoter_motifs/{sample}_{ident}/fimo_cand0.tsv`):
`meme_training.fa` → `discovery/meme_training_dedup.fa` + `discovery/background_order2.txt`
→ `meme/meme.txt` → `cand0.*` → `fimo_cand0.tsv`. Discovery alone: `snakemake --cores 2 motifs_all`.

## Outputs used by `rescue_tss`

`results/promoter_motifs/{sample}_{ident}/`

- **`cand0.bed`**: BED6, the selected `_0` candidates, same columns as the peaks BED
  (window coordinates, name, strand).
- **`fimo_cand0.tsv`**: FIMO `--text` output:
  - `sequence_name` = candidate name, matches column 4 of `cand0.bed`
  - `start` / `stop` = 1-based positions **inside the 55-nt window**
  - windows are oriented 5'→3' along the transcript on both strands
    (window = TSS−50 … TSS+4), so the **TSS is window position 51** on + and −
  - position relative to the TSS = `start − 51` (e.g. a hit at 37–51 spans −14…0)
- `cand0_table.tsv`: every `_0` candidate with TSS, enriched reads, ratio and the
  filter outcome (selected / low_ratio / …), for checking the selection.

FIMO 5.5 reads sequence names that look like coordinates (`X::chr:1-55(+)`, but also
`TSS_POS_4_0` → name `TSS_POS`, start 4) as genomic coordinates and shifts the hit
positions. `scan_candidates` therefore runs FIMO on `cand0_scan.fa` (headers `cand0001`, …)
and maps the IDs back to candidate names with `cand0_ids.tsv` (`scripts/rename_fimo.py`).
`fimo_cand0.tsv` thus always has real candidate names and positions 1–55.
For `scan_background`: check that `sequence_name` in `fimo_random.tsv` matches the
`random.fa` headers and that no `stop` exceeds the window length.

## Config

```yaml
motif background:                 # shared with scan_background
  fimo threshold: 0.001
  fimo options: "--norc"          # windows are already oriented: scan one strand
  # no "motif file": default results/promoter_motifs/{sample}_{ident}/meme/meme.txt

candidate selection:
  phage records: null             # null = {sample}, as in selectAboveThresholdSequences
  min enriched reads: 0
  min enrichment ratio: 1.0
  exclude near _1: 0              # e.g. 55 to drop _0 positions next to a _1 promoter
```

## Open points

- `background_model` / `scan_background` use `config["fasta file"]` (host + phage + spike-in):
  the background composition and ~99% of the random windows then come from the host.
  Restrict both to the phage record(s). Also `2 > {log}` should be `2> {log}`.
- `cand1.bed` (for `rescue_tss`) is not produced yet.
- `scan_background` should read the same motif file:
  `results/promoter_motifs/{sample}_{ident}/meme/meme.txt` (not `test/test_motif.meme`).
- LUZ7 t5 alone gives 5 independent `_1` promoters; MEME-1 is positionally conserved
  but not significant (E = 4.0). Pool `_1` windows over timepoints/replicates next.
