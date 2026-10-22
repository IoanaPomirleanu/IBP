# Transcription unit (TU) module

## Definition

A TU runs from a TSS to a downstream TTS on the same strand, and is kept only if reads start at the TSS and end at the TTS (as in Putzeys et al. 2022, 2023 and the LUZ100 paper).
- One TSS can give several TUs: reads that read through a TTS and stop at a later one.
- A TSS whose reads end at no TTS gets one TU up to its **furthest read 3' end** ("longest read").
- No distance cut-off is needed: the reads decide which TTS belongs to which TSS.

## Rules (`workflow/rules/detectTU.smk`)

| Step | Rule | What it does |
|---|---|---|
| 1 | `tssPoints` | final TSS set (`tss_final.bed`: tier 1 `_1` + tier 2 rescued `_0`), phage only → 1-nt points |
| 2 | `ttsPoints` | TTS windows (`TTS_*_corrected.bed`) → 1-nt points |
| 3 | `readEnds` | primary reads (`samtools view -F 0x904`) of the enriched BAM on the phage → point at the 5' end, 3' end in the name (temporary file) |
| 4 | `readsPerTSS` | `bedtools window -sm -w {TSS window}`: reads starting at each TSS → their 3' ends, labelled with the TSS |
| 5 | `pairTSSandTTS` | `bedtools window -sm -w {TTS window}`: 3' ends near a downstream TTS → reads counted per (TSS, TTS), kept if ≥ `minimum reads` |
| 6 | `fallbackTUs` | TSSs with reads but no pair → TU to the furthest 3' end |
| 7 | `combineTUs` | steps 5 + 6 → final BED + table |
| 8 | `genesPerTU` | genes on the same strand ≥ 90% inside each TU (`bedtools intersect -F 0.9 -s`, GFF) → table + `n_genes`, `genes` |

Steps 1, 2, 6, 7 are awk only; steps 3-5 and 8 use `envs/env_tu.yaml` (bedtools + samtools).

**Phage filter:** steps 1 and 3 keep only contigs whose name contains the sample name (`index($1, "{wildcards.sample}")`, as in `annotateTSS.smk`), so `Sample name` must be part of the phage contig name (e.g. `LUZ7`). Host and spike-in are left out.

**Which TSSs:** tiers 1 and 2 only. Tier 3 (`_0` not rescued) is treated as "not a TSS" and gets no TU. `tss_tiers.tsv` lists all three tiers; `tss_final.bed` is tiers 1 + 2.

**Close TSSs are not merged:** `dedupTraining` (10 nt) only affects the MEME training set. Two TSSs a few nt apart (e.g. LUZ7 100 and 107) and terminal-repeat copies (100 and 74341) each get their own TUs.

## Setup

New files:
- `workflow/rules/detectTU.smk`
- `workflow/envs/env_tu.yaml`

`config/config.yaml` (block in `config/config_TU.yaml`; make sure the line before it ends with a newline):
```
TU detection:
  TSS window: 5        # a read starts "at" a TSS if its 5' end is within this many nt
  TTS window: 20       # a read ends "at" a TTS if its 3' end is within this many nt
  minimum reads: 1     # reads needed to keep a TSS -> TTS pair
  # step 8, genes per TU
  gene annotation: input/genome_data/LUZ7.gff3   # GFF of the phage genes only (e.g. NCBI NC_013691.1)
  gene features: gene                           # GFF column 3 types to use, e.g. "gene" or "gene|ncRNA"
  gene name attribute: Name                     # attribute shown as the gene name, e.g. Name or locus_tag
  gene overlap fraction: 0.9                    # fraction of a gene that must lie inside the TU (papers: 0.9)
```

**GFF:** phage genes only (e.g. the NCBI GFF3 of NC_013691.1 for LUZ7). Its sequence name doesn't have to match the FASTA contig name: step 8 gives every gene the phage contig name of the TUs. Pick `gene name attribute` so the names look like you want (`Name`, `locus_tag`, ...); the paper's `orf01`… names may only match one of them.

`workflow/Snakefile`:
```
include: 'rules/detectTU.smk'
```
and in `rule all`:
```
        expand("results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TU_genes_{sample}_{ident}.tsv", sample=config["Sample name"], ident=config["ID"])
```
(Without a GFF, use `TU_{sample}_{ident}.tsv` instead: steps 1-7 only.)

Step 1 needs `tss_final.bed`, so `rescueTSS` (and the motif rules before it) must be included too; Snakemake runs them first if the file doesn't exist.

## Run

Only the TU table (dry run first):
```
snakemake --use-conda -n results/transcription_units/TU_LUZ7/TU_LUZ7_t5/TU_genes_LUZ7_t5.tsv
snakemake --use-conda --cores 4 results/transcription_units/TU_LUZ7/TU_LUZ7_t5/TU_genes_LUZ7_t5.tsv
```
- After changing a TU rule: add `--forcerun <rule>` (e.g. `--forcerun tssPoints readEnds`); everything after it reruns too.
- If the dry run also lists mapping or peak calling because other rule code changed: add `--rerun-triggers mtime`.
- Runs in seconds to minutes; the BAM is read once.

## Outputs

`results/transcription_units/TU_{sample}/TU_{sample}_{ident}/`

- **`TU_{sample}_{ident}.tsv`**: one line per TU, with header:
  `TU_id, chr, start, end, strand, TSS, TSS_pos, end_id, end_pos, n_reads, type`
  - `start`/`end`: BED coordinates (0-based start); `TSS_pos`/`end_pos`: 1-based positions
  - `end_id`: the TTS name, or `longest_read`
  - `n_reads`: reads from the TSS ending at that TTS (`TSS_TTS`), or all reads from the TSS (`longest_read`)
  - `type`: `TSS_TTS` or `longest_read`: why the TU is there
- **`TU_genes_{sample}_{ident}.tsv`** (step 8): the same table plus `n_genes` and `genes` (names joined with ` | `, `-` if none), like Table S4 of the 2022 paper
- **`TU_{sample}_{ident}.bed`**: same TUs as BED6 for IGV (score = `n_reads`)
- intermediate files for checking: `TSS_points_*`, `TTS_points_*`, `TSS_read_ends_*`, `TU_pairs_*`, `TU_longest_read_*`

## Checks

```
cd results/transcription_units/TU_LUZ7/TU_LUZ7_t5
cut -f2 TU_LUZ7_t5.tsv | sort | uniq -c        # only LUZ7 lines
cut -f11 TU_LUZ7_t5.tsv | sort | uniq -c       # TSS_TTS vs longest_read
awk '$1=="LUZ7"' TTS_points_LUZ7_t5.bed        # which TTSs were called
```
Pipeline TSS = published TSS + 1.

## Result on LUZ7 t5 (tier 1 only, 12 Oct 2026)

Compared with Supplementary Table S4 of Putzeys et al. 2022 (86 TUs, t5 + t10 + t20 pooled):
- 17 LUZ7 TUs; 14 distinct after removing the TSS at 107 (shoulder of PE1) and the terminal-repeat copies. **All 14 are in Table S4**, ends within 4 nt.
- The 6 early promoters found at t5 (PE1, 3, 4, 5, 7, 8) have 22 TUs in the paper: 14 found. The 8 missing ones end at TTSs not called at t5 (1924, 2649, 4350, 9964, 10511).
- The other 64 paper TUs start at promoters not found at t5 (PE2, PE6, PE9 and the middle/late promoters).

Not yet done: rerun with tiers 1 + 2 (`tss_final.bed`) and genes per TU (step 8, needs the LUZ7 GFF), t10/t20.
