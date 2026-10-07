# rescue_tss
# Decide which _0 TSS candidates to rescue based on promoter motifs.
# Motif p-value cutoffs are calibrated on random genomic windows (scan_background),
# so the false-positive rate of the rescue is measured, not assumed.
# Output: final table with a tier column (1 = _1 TSS, 2 = rescued _0, 3 = _0 not rescued).


rule rescue_tss:
    input:
        random_fimo = "results/promoter_motifs/{sample}_{ident}/background/fimo_random.tsv",
        random_fa   = "results/promoter_motifs/{sample}_{ident}/background/random.fa",
        cand0_fimo  = "results/promoter_motifs/{sample}_{ident}/fimo_cand0.tsv",  # <-- CHANGE: FIMO on _0 (Zinaida, scan_candidates)
        cand0_bed   = "results/promoter_motifs/{sample}_{ident}/cand0.bed",       # <-- CHANGE: _0 candidates (Ioana)
        cand1_bed   = "results/promoter_motifs/{sample}_{ident}/cand1.bed"        # <-- CHANGE: _1 TSS (Ioana)
    output:
        tiers   = "results/promoter_motifs/{sample}_{ident}/rescue/tss_tiers.tsv",
        bed     = "results/promoter_motifs/{sample}_{ident}/rescue/tss_final.bed",
        rates   = "results/promoter_motifs/{sample}_{ident}/rescue/background_rates.tsv",
        cutoffs = "results/promoter_motifs/{sample}_{ident}/rescue/calibrated_cutoffs.tsv",
        summary = "results/promoter_motifs/{sample}_{ident}/rescue/rescue_summary.txt"
    params:
        script     = workflow.basedir + "/scripts/rescue_tss.py",
        target_fpr = config["motif rescue"]["target fpr"],
        position   = rescue_position_args,
        outdir     = "results/promoter_motifs/{sample}_{ident}/rescue"
    log:
        "logs/rescue_tss/{sample}_{ident}.log"
    conda:
        "../envs/env_scan_background.yaml"
    shell:
        """
        python {params.script} \
            --random-fimo {input.random_fimo} --random-fa {input.random_fa} \
            --cand0-fimo {input.cand0_fimo} \
            --cand0-bed {input.cand0_bed} --cand1-bed {input.cand1_bed} \
            --target-fpr {params.target_fpr} {params.position} \
            --outdir {params.outdir} > {log} 2>&1
        """
