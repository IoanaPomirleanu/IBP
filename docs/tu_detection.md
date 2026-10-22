# Transcription units (TUs)

Rules: `workflow/rules/detectTU.smk`, environment `workflow/envs/env_tu.yaml` (bedtools + samtools).

## Definition

A TU runs from a TSS to a downstream TTS on the same strand, and is kept only if reads start at the TSS
and end at the TTS (as in the LUZ7, LUZ100 and six-phage ONT-cappable-seq papers).

- One TSS can give several TUs: reads that read through a TTS and stop at a later one.
- A TSS whose reads end at no TTS gets one TU up to its **furthest read 3' end** ("longest read").
- No distance cut-off is needed: the reads decide which TTS belongs to which TSS.
- Genes belong to a TU when at least 90% of the gene lies inside it, on the same strand (as in the papers,
  `bedtools intersect -F 0.9 -s`).

## Rules

| Step | Rule | What it does |
|---|---|---|
| 1 | `tssPoints` | final TSS set (`tss_final.bed`: tiers 1 + 2), phage contig only → 1-nt points |
| 2 | `ttsPoints` | TTS windows (`TTS_*_corrected.bed`) → 1-nt points |
| 3 | `readEnds` | primary reads of the enriched BAM on the phage (`samtools view -F 0x904`) → a point at the 5' end, with the 3' end kept (temporary file) |
| 4 | `readsPerTSS` | reads whose 5' end is within `TSS window` nt of a TSS → their 3' ends, labelled with that TSS |
| 5 | `pairTSSandTTS` | 3' ends within `TTS window` nt of a downstream TTS → reads counted per (TSS, TTS); pairs with ≥ `minimum reads` are TUs |
| 6 | `fallbackTUs` | TSSs with reads but no pair → one TU to the furthest 3' end |
| 7 | `combineTUs` | steps 5 + 6 → table and BED |
| 8 | `genesPerTU` | genes on the same strand ≥ `gene overlap fraction` inside each TU (GFF3) → table with `n_genes`, `genes` |

Steps 1, 2, 6 and 7 are awk only.

**Phage filter:** steps 1 and 3 keep only contigs whose name contains the `Sample name`, so host and
spike-in are left out.

**Which TSSs:** tiers 1 and 2 only. Tier 3 (`_0` not rescued) gets no TU.

**Close TSSs are not merged:** two TSSs a few nt apart and terminal-repeat copies each get their own TUs
(`dedup distance` only affects motif training).

**GFF:** phage genes only. Its sequence name does not have to match the FASTA: step 8 gives every gene the
phage contig name of the TUs. Choose `gene name attribute` (`Name`, `locus_tag`, ...) for the names you want.

## Config

The `TU detection` block in `config/config.yaml`: `TSS window`, `TTS window`, `minimum reads`,
`gene annotation`, `gene features`, `gene name attribute`, `gene overlap fraction`.

## Outputs

`results/transcription_units/TU_{sample}/TU_{sample}_{ID}/`

- **`TU_genes_{sample}_{ID}.tsv`**: one line per TU:
  `TU_id, chr, start, end, strand, TSS, TSS_pos, end_id, end_pos, n_reads, type, n_genes, genes`
  - `start` / `end`: BED coordinates (0-based start); `TSS_pos` / `end_pos`: 1-based positions
  - `end_id`: the TTS name, or `longest_read`
  - `n_reads`: reads from the TSS ending at that TTS (`TSS_TTS`), or all reads from the TSS (`longest_read`)
  - `type`: `TSS_TTS` or `longest_read`
  - `genes`: gene names joined with ` | `, `-` if none
- **`TU_{sample}_{ID}.tsv`**: the same without genes (steps 1-7)
- **`TU_{sample}_{ID}.bed`**: BED6 for IGV (score = `n_reads`)
- intermediate files for checking: `TSS_points_*`, `TTS_points_*`, `TSS_read_ends_*`, `TU_pairs_*`, `TU_longest_read_*`

## Run only this module

```bash
snakemake --use-conda -n --configfile config/<dataset>.yaml \
    results/transcription_units/TU_<sample>/TU_<sample>_<ID>/TU_genes_<sample>_<ID>.tsv
```

- After changing a TU rule: add `--forcerun <rule>` (everything after it reruns too).
- If the dry run also lists mapping or peak calling because other rule code changed:
  add `--rerun-triggers mtime`.

## Checks

```bash
cd results/transcription_units/TU_<sample>/TU_<sample>_<ID>
cut -f2 TU_<sample>_<ID>.tsv | sort | uniq -c      # only phage lines
cut -f11 TU_<sample>_<ID>.tsv | sort | uniq -c     # TSS_TTS vs longest_read
```

Many `longest_read` TUs usually mean TTSs that were not called; check the TTS settings.
