# rescue_tss
# Decide which _0 TSS candidates to rescue based on promoter motifs.
# Motif p-value cutoffs are calibrated on random genomic windows (scan_background),
# so the false-positive rate of the rescue is measured, not assumed.
# Output: final table with a tier column (1 = _1 TSS, 2 = rescued _0, 3 = _0 not rescued).


rule rescue_tss:
    input:
        random_fimo = "results/promoter_motifs/{sample}_{ident}/background/fimo_random.tsv",
        random_fa   = "results/promoter_motifs/{sample}_{ident}/background/random.fa",
        cand0_fimo  = rules.scan_candidates.output.fimo,  # takes whatever file scan_candidates writes as it's fimo output
        cand0_bed   = rules.select_cand0.output.bed,
        tss_fa = SC_TSS_DIR + "/TSS_seq_{sample}_enriched_{ident}.fa.out",
        
    output:
        tss_tiers   = "results/promoter_motifs/{sample}_{ident}/rescue/tss_tiers.tsv",
        bed     = "results/promoter_motifs/{sample}_{ident}/rescue/tss_final.bed",
        rates   = "results/promoter_motifs/{sample}_{ident}/rescue/background_rates.tsv",
        cutoffs = "results/promoter_motifs/{sample}_{ident}/rescue/calibrated_cutoffs.tsv",
        summary = "results/promoter_motifs/{sample}_{ident}/rescue/rescue_summary.txt"
    params:
        script     = workflow.basedir + "/scripts/rescue_tss.py",
        target_fpr = config["motif rescue"]["target fpr"],
        outdir     = "results/promoter_motifs/{sample}_{ident}/rescue",
        spike = config["spike-in name"]
    log:
        "logs/rescue_tss/{sample}_{ident}.log"
    conda:
        "../envs/env_scan_background.yaml"
    shell:
        """
        python {params.script} \
            --random-fimo {input.random_fimo} --random-fa {input.random_fa} \
            --cand0-fimo {input.cand0_fimo} \
            --cand0-bed {input.cand0_bed} --tss-fa {input.tss_fa}\
            --target-fpr {params.target_fpr} \
            --outdir {params.outdir} > {log} 2>&1
        """
