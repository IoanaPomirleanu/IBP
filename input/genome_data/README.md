# Reference files

Files here are not tracked by git. Put their paths in your config (`fasta file`, `gene annotation`).

**FASTA** (`fasta file`): the phage, host and spike-in sequences concatenated in one file, e.g.
`PHAGE_host_spikein.fasta`:

```bash
cat PHAGE.fasta HOST.fasta SPIKEIN.fasta > PHAGE_host_spikein.fasta
```

- The phage header must contain the `Sample name` (e.g. `>LUZ7`): the motif and TU steps use it to
  keep only the phage.
- The spike-in header must contain the `spike-in name` (e.g. `>FLuc`): its TSS sets the threshold.

**GFF3** (`gene annotation`, e.g. `PHAGE.gff3`): the phage genes only, for the genes per TU. The
sequence name in the GFF does not have to match the FASTA header.
