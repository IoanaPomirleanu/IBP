# ONT-cappable-seq pipeline extensions

Improvements to the ONT-cappable-seq workflow (Putzeys et al., 2022) for mapping
transcription start/termination sites in Pseudomonas phages from nanopore
long-read data.

## Goals
- **Dynamic TSS calling**: replace the fixed, spike-in-derived enrichment-ratio
  threshold with a data-driven decision that integrates promoter motif evidence
  (σ70-like and phage-specific motifs).
- **Automated TU module**: build transcription units from TSS/TTS pairs and
  read-spanning evidence, including multiple TUs per TSS (readthrough), replacing manual curation.

## Status
Work in progress.

## Based on
This project builds on the [ONT-cappable-seq pipeline](https://github.com/LoGT-KULeuven/ONT-cappable-seq) by LoGT-KULeuven, starting from revision [`3920110`](https://github.com/LoGT-KULeuven/ONT-cappable-seq/tree/392011099533cfb797612c080c6ea5c8693e0c7d). See the original repository for installation and usage instructions.