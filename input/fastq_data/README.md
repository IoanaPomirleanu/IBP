# Sequencing reads

Files here are not tracked by git. Put their paths in your config (`enriched fastq`, `control fastq`).

One enriched and one control library per sample, e.g.:

```
PHAGE_enriched.fastq
PHAGE_control.fastq
```

For several time points or replicates, use one config per sample with its own `ID` (no `_` in
`Sample name` or `ID`).
