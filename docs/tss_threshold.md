# TSS candidates and the spike-in threshold

Rules: `workflow/rules/annotateTSS.smk` (after the peak calling of the original pipeline).

## Idea

In ONT-cappable-seq, primary transcripts (5' triphosphate ends, captured by Cappable-seq) are
enriched in the enriched library relative to the control library. A TSS is called when its
**enrichment ratio** (enriched / control reads at the peak) is high enough.

The original pipeline used a fixed `TSS Threshold` from the config, taken from the paper per dataset.
Here the threshold is read from the data itself: an in vitro transcribed **spike-in RNA** with a known
5' triphosphate end is added to each sample, so its 5' end behaves like a true primary TSS. Its enrichment ratio is
used as the threshold for the phage TSSs of that sample.

Nothing is thrown away: every candidate stays in the output with a label, so that weaker candidates can
still be rescued later by their promoter motif ([motifs.md](motifs.md)).

## Rules

| Rule | What it does |
|---|---|
| `getPeaks`, `splitOverlappingPeaks`, `getOverlappingPeaksWithError`, `getAllOverlappingPeaks` | pair each enriched 5' peak with a control peak at the same position, or within `peak alignment error` nt |
| `getInformationFromOverlappingPeaks` | one row per paired peak with both counts and the enrichment ratio (`enr_ratios_*.plus.csv` / `.minus.csv`); no threshold filter |
| `findSpikeInThreshold` | finds the + strand candidate on the contig whose name contains `spike-in name`; if several, the one with most enriched reads. Its ratio is the threshold (`spike_in_threshold_*.txt`: contig, position, ratio). Stops with "Spike-in not found in data" if there is none |
| `posPeaksToBedFile`, `negPeaksToBedFile` | promoter window per candidate (`upstream` / `downstream` nt around the TSS); name `TSS_POS_<n>_1` / `TSS_NEG_<n>_1` if ratio >= threshold, `_0` if below |
| `combinePosAndNegBedFiles` | both strands in one BED; window starts below 0 (TSS near the genome start) are set to 0 |
| `extractSequencesTSS` | window sequences (`TSS_seq_*.fa.out`), oriented 5'->3' on both strands |

## Config

```yaml
spike-in name: FLuc        # part of the spike-in contig name
peak alignment error: 1    # nt allowed between an enriched and a control peak
minimum coverage:          # reads needed at a peak (enriched / control)
TSS sequence extraction:   # promoter window: upstream / downstream of the TSS
```

## Outputs

`results/transcript_boundaries/TSS_{sample}/TSS_{sample}_{ID}/`

- `enr_ratios_{sample}_{ID}.plus.csv` / `.minus.csv`: all paired candidates with counts and ratio
- `spike_in_threshold_{sample}_{ID}.txt`: the threshold
- `{sample}_enriched_{ID}_peaks.bed`: all candidates with their `_1` / `_0` label
- `TSS_seq_{sample}_enriched_{ID}.fa.out`: their promoter windows (input for the motif module)

## Checks

- `spike_in_threshold_*.txt` exists and names the spike-in contig. If the run stops with
  "Spike-in not found", check `spike-in name` against the FASTA headers and that the spike-in has reads.
- Count the labels on the phage contig: `awk '$1 ~ /PHAGE/' *_enriched_*_peaks.bed | grep -c '_1$'`.
- TSS position in a window: `BED end - downstream` (+ strand), `BED start + downstream` (- strand).

## Limitations

- A TSS needs a peak in the control library too. Sites with reads only in the enriched library are not
  candidates. Lowering `minimum coverage: control: TSS` or raising `peak alignment error` recovers some.
- Reads in terminal repeats map to both copies (MAPQ 0), so the ratio of each copy is lower.
