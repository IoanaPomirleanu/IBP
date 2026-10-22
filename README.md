# ONT-cappable-seq for phage transcriptomics

Extensions to the [ONT-cappable-seq pipeline](https://github.com/LoGT-KULeuven/ONT-cappable-seq)
(Putzeys et al., 2022) for mapping the transcriptional landscape of phages from nanopore
long-read data: transcription start sites (TSS), termination sites (TTS) and transcription units (TU).

## What this fork adds

- **TSS threshold from the data.** The enrichment ratio of the spike-in RNA's own TSS is used as the
  threshold, instead of a fixed value in the config. All candidates are kept and labelled
  `_1` (above the threshold) or `_0` (below).
- **Promoter motifs.** MEME finds motifs in the promoter windows of the `_1` TSSs. FIMO then scans the
  `_0` candidates, and a `_0` with a motif at the right distance from its TSS is rescued. The p-value
  cut-offs are calibrated on random genome windows, so the false-positive rate of the rescue is measured.
- **Known motifs.** Motifs from the literature (a MEME file or IUPAC consensus sequences) can be added
  to the rescue, alone or together with the de novo motifs.
- **Transcription units.** TSS and TTS are paired by the reads that span them; a TSS without a TTS gets
  a TU up to its longest read. The genes inside each TU are listed.

## Pipeline overview

| Step | Rules | Explained in |
|---|---|---|
| 1. Read processing and mapping | `mapReads.smk` | original pipeline |
| 2. 5' and 3' peak calling | `determineTranscriptionBoundaries.smk` | original pipeline |
| 3. TSS candidates and spike-in threshold | `annotateTSS.smk` | [docs/tss_threshold.md](docs/tss_threshold.md) |
| 4. TTS | `annotateTTS.smk` | original pipeline |
| 5. Motif discovery, candidate scan, rescue, known motifs | `motifs.smk`, `scanCandidates.smk`, `scanBackground.smk`, `knownMotifs.smk`, `rescueTSS.smk` | [docs/motifs.md](docs/motifs.md) |
| 6. Transcription units and genes | `detectTU.smk` | [docs/tu_detection.md](docs/tu_detection.md) |

## Input

In `input/genome_data/` and `input/fastq_data/` (not tracked by git):

- **FASTA**: phage + host + spike-in sequences concatenated in one file. The phage contig name must
  contain the `Sample name` (e.g. `>LUZ7`), and the spike-in contig name must contain `spike-in name`
  (e.g. `FLuc`).
- **FASTQ**: the enriched and the control library of one sample.
- **GFF3** of the phage genes only (e.g. from NCBI), for the genes per TU.

## How to run

The pipeline uses Snakemake with conda environments (`workflow/envs/`).

1. `config/config.yaml` holds the defaults for every setting, with a short comment per key.
2. For each dataset, write a small config with only the keys that differ: `Sample name`, `ID`, the
   input paths, and any changed settings. Use no `_` in `Sample name` or `ID`.
3. Run a dry run first, then the pipeline:

```bash
snakemake --use-conda -n --configfile config/<dataset>.yaml
snakemake --use-conda --cores 8 --configfile config/<dataset>.yaml
```

## Main outputs

Paths per `{sample}_{ID}`:

| File | Content |
|---|---|
| `results/transcript_boundaries/TSS_{sample}/.../spike_in_threshold_*.txt` | spike-in contig, TSS position and enrichment ratio (= threshold) |
| `results/promoter_motifs/{sample}_{ID}/meme/meme.txt` | de novo motifs |
| `results/promoter_motifs/{sample}_{ID}/rescue/tss_tiers.tsv` | every TSS candidate with its tier: 1 = `_1`, 2 = rescued `_0`, 3 = `_0` not rescued |
| `results/promoter_motifs/{sample}_{ID}/rescue/tss_final.bed` | final TSS set (tiers 1 + 2) |
| `results/promoter_motifs/{sample}_{ID}/rescue/rescue_summary.txt` | motifs used, cut-offs, expected false rescues |
| `results/transcript_boundaries/TTS_{sample}/.../TTS_*_corrected.bed` | TTS windows |
| `results/transcription_units/TU_{sample}/.../TU_genes_*.tsv` | TUs with their genes |
| `results/transcription_units/TU_{sample}/.../TU_*.bed` | TUs for IGV |

## Credits

Based on [LoGT-KULeuven/ONT-cappable-seq](https://github.com/LoGT-KULeuven/ONT-cappable-seq) at
revision [`3920110`](https://github.com/LoGT-KULeuven/ONT-cappable-seq/tree/392011099533cfb797612c080c6ea5c8693e0c7d).
Method and benchmark data:

- Putzeys L, Boon M, Lammens E-M, Kuznedelov K, Severinov K, Lavigne R. Development of ONT-cappable-seq to
  unravel the transcriptional landscape of *Pseudomonas* phages. *Comput Struct Biotechnol J* 2022;20:2624–2638.
  [doi:10.1016/j.csbj.2022.05.034](https://doi.org/10.1016/j.csbj.2022.05.034) (LUZ7)
- Putzeys L, Poppeliers J, Boon M, Lood C, Vallino M, Lavigne R. Transcriptomics-driven characterization of
  LUZ100, a T7-like *Pseudomonas* phage with temperate features. *mSystems* 2023;8(2):e01189-22.
  [doi:10.1128/msystems.01189-22](https://doi.org/10.1128/msystems.01189-22) (LUZ100)
- Putzeys L, Wicke L, Boon M, van Noort V, Vogel J, Lavigne R. Refining the transcriptional landscapes for
  distinct clades of virulent phages infecting *Pseudomonas aeruginosa*. *microLife* 2024;5:uqae002.
  [doi:10.1093/femsml/uqae002](https://doi.org/10.1093/femsml/uqae002) (14-1, LUZ19, LUZ24, PAK_P3, YuA, phiKZ)

Developed by: Zinaida Barseghyan, Lisa Gillebert and Ioana Pomirleanu.
