# Promoter motifs: discovery, candidate scan and rescue

Rules: `motifs.smk`, `scanBackground.smk`, `scanCandidates.smk`, `knownMotifs.smk`,
`rescueTSS.smk`, plus `selectAboveThresholdSequences` at the end of `annotateTSS.smk`. Scripts in `workflow/scripts/`.

## Idea

The spike-in threshold ([tss_threshold.md](tss_threshold.md)) splits the TSS candidates into `_1` (above)
and `_0` (below). Some `_0` are real promoters with weak expression or a weak control peak. The `_1`
promoters show which sequence motifs this phage uses, so a `_0` that carries such a motif **at the same
distance from its TSS** is likely real too and is rescued.

To keep the rescue honest, the motif p-value cut-offs are calibrated on random windows of the phage
genome: the cut-offs are chosen so that at most `target fpr` (e.g. 5%) of random windows would be
rescued. The false-positive rate is therefore measured, not assumed.

All windows (`_1`, `_0`, random) have the same length and are oriented 5'->3' along the transcript, so a
motif position in the window means the same distance to the TSS everywhere.

## Steps and rules

| Step | Rule | What it does |
|---|---|---|
| 1. Training set | `selectAboveThresholdSequences` | promoter windows of the `_1` TSSs on the phage contig → `meme_training.fa` |
| | `dedupTraining` | `_1` windows less than `dedup distance` nt apart, or identical (terminal-repeat copies), count as one promoter (training only; the TSS calls are not changed) |
| 2. Discovery | `motifBackgroundModel` | Markov background (order `markov order`) of the phage genome only |
| | `discoverMotifs` | MEME (`zoops`, given strand only) on the deduplicated windows → `meme/meme.txt`; skipped if fewer than `min sequences` windows |
| 3. Known motifs (optional) | `combineMotifs` | de novo motifs + known motifs (named `KNOWN-<name>`) in one file, `motifs/scan_motifs.meme` |
| | `scanTraining` | FIMO on the `_1` windows: where each known motif sits relative to the TSS |
| 4. Background | `backgroundModel` | base composition of the phage genome (`bg.txt`, for FIMO) |
| | `scanBackground` | `n random windows` random phage windows (same length, random strand), FIMO with the same settings → `background/fimo_random.tsv` |
| 5. Candidates | `selectCand0` | the `_0` candidates on the phage contig, optional filters (`min enriched reads`, `min enrichment ratio`, `exclude near _1`) → `cand0.bed`, `cand0_table.tsv` |
| | `scanCandidates` | FIMO on the `_0` windows, same motifs and settings → `fimo_cand0.tsv` |
| 6. Rescue | `rescueTSS` | motif selection, position rule, calibration, rescue → `rescue/` |

## Rescue in detail (`rescueTSS`)

1. **Motif selection.** A motif is used when its MEME sites in the `_1` windows are positionally conserved:
   at least `min training sites` sites, whose start positions lie within `max site spread` nt of each
   other. A real promoter element sits at a fixed distance from the TSS; a motif scattered over the window
   is probably not a promoter element. Optional: `max evalue` (MEME E-value filter) or `motifs`
   (a fixed list, e.g. `[MEME-1]`).
2. **Position rule.** A FIMO hit only counts if it starts where the motif was found in the `_1` windows,
   ± `position slack` nt (`position rule: auto`). `off` = anywhere in the window; `[26, 41]` = fixed range.
   The same rule is applied to the random windows.
3. **Calibration.** Starting at `fimo threshold`, each motif's p-value cut-off is lowered until at most
   `target fpr` of the random windows have a hit; then all cut-offs are lowered together until at most
   `target fpr` of the random windows have a hit of **any** used motif. `target fpr` is therefore the
   false-positive rate of the whole rescue.
4. **Rescue.** A `_0` candidate is rescued if any used motif hits it inside its position range with a
   p-value at or below that motif's cut-off.
5. **Tiers.** 1 = `_1` (above the spike-in threshold), 2 = rescued `_0`, 3 = `_0` not rescued.

**Known motifs** are always used. Their position comes from `positions:` in the config (motif start
relative to the TSS, TSS = 0), or is learned from their hits in the `_1` windows (same site-count and
spread rules). If neither works, they count anywhere in the window, which rarely rescues anything, so
give a position. With `use de novo motifs: false`, MEME is not run and only known motifs are used (useful
when there are too few `_1` promoters).

## Config

The blocks `motifs`, `motif background`, `candidate selection`, `motif rescue` and `known motifs` in
`config/config.yaml`, each key with a short comment.

## Outputs

`results/promoter_motifs/{sample}_{ID}/`

| File | Content |
|---|---|
| `meme_training.fa`, `discovery/meme_training_dedup.fa` | `_1` windows, before and after deduplication |
| `meme/meme.txt` | de novo motifs with E-values and training sites |
| `motifs/scan_motifs.meme`, `motifs/known_motifs.tsv` | combined motif file and list of known motifs (only with known motifs) |
| `background/fimo_random.tsv` | FIMO hits in the random windows |
| `cand0_table.tsv` | every `_0` candidate and why it was selected or not (`selected`, `low_ratio`, `near_1`, ...) |
| `fimo_cand0.tsv` | FIMO hits in the `_0` windows (positions 1..window length) |
| `rescue/motif_selection.tsv` | per motif: source, E-value, training-site starts, allowed position range, used or not and why |
| `rescue/calibrated_cutoffs.tsv`, `rescue/background_rates.tsv` | the calibration steps |
| `rescue/tss_tiers.tsv` | every candidate with tier, best motif and p-value |
| `rescue/tss_final.bed` | final TSS set: tiers 1 + 2 (input of the TU module) |
| `rescue/rescue_summary.txt` | motifs used, background rate, `_0` scanned and rescued, expected false rescues |

## Checks

- `rescue/motif_selection.tsv`: are the motifs you expect (e.g. the σ70 -10 box) used, at a sensible
  distance from the TSS (`allowed_start_rel_tss`)?
- `rescue/rescue_summary.txt`: expected false rescues compared with the number rescued.
- `cand0_table.tsv`: whether a known promoter was a `_0` candidate at all, and if not, why.

## Notes

- FIMO reads sequence names that look like coordinates (e.g. `TSS_POS_4_0`) as positions and shifts the
  hits, so the scans use short IDs (`cand0001`, `trn0001`) that are mapped back to candidate names.
- On the + strand the TSS is window position `upstream + 1`, on the - strand position `upstream`
  (1-nt shift from the window construction). The position slack absorbs it.
- With few `_1` promoters (e.g. one early time point) MEME has little to learn from; motifs may be
  non-significant or fail the positional test. Pooling time points or using known motifs helps.
