# ONT-cappable-seq pipeline extensions

Improvements to the ONT-cappable-seq workflow (Putzeys et al., 2022) for mapping
transcription start/termination sites in Pseudomonas phages from nanopore
long-read data.

## Goals
- **Dynamic TSS calling**: replace the fixed, spike-in-derived enrichment-ratio
  threshold with a data-driven decision that integrates promoter motif evidence
  (σ70-like and phage-specific motifs).
- **Automated TU module**: build transcription units from TSS/TTS pairs and
  read-spanning evidence, including multiple TUs per TSS (readthrough) and
  operon-overlap detection, replacing manual curation.

## Status
Work in progress.
