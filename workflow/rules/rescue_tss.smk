# rescue_tss
# Decide which _0 TSS candidates to rescue based on promoter motifs.
# Only positionally conserved motifs are used, and a hit only counts at the
# position where MEME found the motif in the _1 promoters (motif_selection.tsv).
# Motif p-value cutoffs are calibrated on random genomic windows (scan_background)
# with the same position rule, so the false-positive rate of the rescue as a
# whole (any motif) is measured, not assumed.
# Output: final table with a tier column (1 = _1 TSS, 2 = rescued _0, 3 = _0 not rescued).

RESCUE = config.get("motif rescue", {})


def _rescue_opts(wc):
    o = []
    if RESCUE.get("motifs"):
        o.append("--motifs " + " ".join(RESCUE["motifs"]))
    if RESCUE.get("max evalue") is not None:
        o.append(f"--max-evalue {RESCUE['max evalue']}")
    pos = RESCUE.get("position rule", "auto")
    if pos is None or pos is False or pos == "off":
        o.append("--position off")
    elif isinstance(pos, (list, tuple)):
        o.append(f"--position {pos[0]} {pos[1]}")
    else:
        o.append("--position auto")
    o.append(f"--position-slack {RESCUE.get('position slack', 3)}")
    o.append(f"--min-training-sites {RESCUE.get('min training sites', 3)}")
    o.append(f"--max-site-spread {RESCUE.get('max site spread', 12)}")
    return " ".join(o)


rule rescue_tss:
    input:
        random_fimo = "results/promoter_motifs/{sample}_{ident}/background/fimo_random.tsv",
        random_fa   = "results/promoter_motifs/{sample}_{ident}/background/random.fa",
        cand0_fimo  = rules.scan_candidates.output.fimo,  # takes whatever file scan_candidates writes as it's fimo output
        cand0_bed   = rules.select_cand0.output.bed,
        tss_fa      = SC_TSS_DIR + "/TSS_seq_{sample}_enriched_{ident}.fa.out",
        motifs      = SC_MOTIF_FILE,                       # training sites + E-values
    output:
        tss_tiers = "results/promoter_motifs/{sample}_{ident}/rescue/tss_tiers.tsv",
        bed       = "results/promoter_motifs/{sample}_{ident}/rescue/tss_final.bed",
        rates     = "results/promoter_motifs/{sample}_{ident}/rescue/background_rates.tsv",
        cutoffs   = "results/promoter_motifs/{sample}_{ident}/rescue/calibrated_cutoffs.tsv",
        selection = "results/promoter_motifs/{sample}_{ident}/rescue/motif_selection.tsv",
        summary   = "results/promoter_motifs/{sample}_{ident}/rescue/rescue_summary.txt"
    params:
        script     = workflow.basedir + "/scripts/rescue_tss.py",
        target_fpr = RESCUE.get("target fpr", 0.05),
        start_cut  = config["motif background"]["fimo threshold"],   # = FIMO --thresh of both scans
        up         = config["TSS sequence extraction"]["upstream"],
        down       = config["TSS sequence extraction"]["downstream"],
        opts       = _rescue_opts,
        outdir     = "results/promoter_motifs/{sample}_{ident}/rescue",
    log:
        "logs/rescue_tss/{sample}_{ident}.log"
    conda:
        "../envs/env_scan_background.yaml"
    shell:
        """
        python {params.script} \
            --random-fimo {input.random_fimo} --random-fa {input.random_fa} \
            --cand0-fimo {input.cand0_fimo} \
            --cand0-bed {input.cand0_bed} --tss-fa {input.tss_fa} \
            --motif-file {input.motifs} \
            --target-fpr {params.target_fpr} --start-cutoff {params.start_cut} \
            --up {params.up} --down {params.down} {params.opts} \
            --outdir {params.outdir} > {log} 2>&1
        """
