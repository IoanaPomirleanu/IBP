# =============================================================================
# select_cand0 + scan_candidates
#
# select_cand0
#   Pick the _0 TSS candidates (peaks below the spike-in threshold) on the
#   phage genome. Their windows are taken from the same TSS_seq_*.fa.out file
#   as meme_training.fa, so _0 and _1 windows have identical geometry.
#   Optional floors on enriched reads / enrichment ratio limit how many _0
#   candidates are scanned: expected false rescues = background hit rate x
#   number of _0 candidates scanned.
#
# scan_candidates
#   FIMO on the _0 windows with exactly the same motif file, background file,
#   threshold and options as scan_background, so the cut-offs calibrated on
#   the random windows apply directly to these hits (rescue_tss).
# =============================================================================

SC_MB = config["motif background"]
SC_CS = config.get("candidate selection", {})

SC_TSS_DIR = "results/transcript_boundaries/TSS_{sample}/TSS_{sample}_{ident}"
SC_PM_DIR = "results/promoter_motifs/{sample}_{ident}"

# one motif file for scan_background and scan_candidates; may contain
# {sample} / {ident}. Default: the MEME output of discover_motifs.
SC_MOTIF_FILE = SC_MB.get("motif file") or SC_PM_DIR + "/meme/meme.txt"


def _phage_records(wc):
    recs = SC_CS.get("phage records")
    return " ".join(recs) if recs else wc.sample


rule select_cand0:
    input:
        peaks=SC_TSS_DIR + "/{sample}_enriched_{ident}_peaks.bed",
        seqs=SC_TSS_DIR + "/TSS_seq_{sample}_enriched_{ident}.fa.out",
        ratios=[SC_TSS_DIR + "/enr_ratios_{sample}_{ident}.plus.csv",
                SC_TSS_DIR + "/enr_ratios_{sample}_{ident}.minus.csv"],
    output:
        bed=SC_PM_DIR + "/cand0.bed",
        fa=SC_PM_DIR + "/cand0.fa",
        scan_fa=SC_PM_DIR + "/cand0_scan.fa",
        ids=SC_PM_DIR + "/cand0_ids.tsv",
        table=SC_PM_DIR + "/cand0_table.tsv",
    params:
        script=workflow.basedir + "/scripts/select_cand0.py",
        records=_phage_records,
        min_reads=SC_CS.get("min enriched reads", 0),
        min_ratio=SC_CS.get("min enrichment ratio", 1.0),
        near1=SC_CS.get("exclude near _1", 0),
    log:
        "logs/select_cand0/{sample}_{ident}.log",
    conda:
        "../envs/env_scan_background.yaml"
    shell:
        "python {params.script} --peaks {input.peaks} --seqs {input.seqs} "
        "--ratios {input.ratios} --records {params.records} "
        "--min-reads {params.min_reads} --min-ratio {params.min_ratio} "
        "--exclude-near-1 {params.near1} "
        "--out-bed {output.bed} --out-fa {output.fa} --out-scan-fa {output.scan_fa} "
        "--out-ids {output.ids} --out-table {output.table} 2> {log}"


rule scan_candidates:
    input:
        fa=rules.select_cand0.output.scan_fa,   # short IDs: FIMO can't misread them
        ids=rules.select_cand0.output.ids,
        motifs=SC_MOTIF_FILE,
        bg=SC_PM_DIR + "/bg.txt",                 # from background_model
    output:
        fimo=SC_PM_DIR + "/fimo_cand0.tsv",
    params:
        rename=workflow.basedir + "/scripts/rename_fimo.py",
        raw=SC_PM_DIR + "/fimo_cand0.raw.tsv",
        thresh=SC_MB["fimo threshold"],
        fimo_opts=SC_MB["fimo options"],
    log:
        "logs/scan_candidates/{sample}_{ident}.log",
    conda:
        "../envs/env_scan_background.yaml"
    shell:
        """
        (
        if grep -q '^MOTIF' {input.motifs} && [ -s {input.fa} ]; then
            echo "scanning $(grep -c '>' {input.fa}) _0 windows with $(grep -c '^MOTIF' {input.motifs}) motif(s)" >&2
            # identical settings to scan_background
            fimo --text {params.fimo_opts} --thresh {params.thresh} --bgfile {input.bg} \
                {input.motifs} {input.fa} > {params.raw}
            # window IDs -> candidate names (positions stay 1..window length)
            python {params.rename} {params.raw} {input.ids} {output.fimo}
            rm -f {params.raw}
        else
            echo "no motifs in {input.motifs} or no _0 windows in {input.fa}: writing an empty scan" >&2
            printf 'motif_id\\tmotif_alt_id\\tsequence_name\\tstart\\tstop\\tstrand\\tscore\\tp-value\\tq-value\\tmatched_sequence\\n' > {output.fimo}
        fi
        ) 2> {log}
        """
