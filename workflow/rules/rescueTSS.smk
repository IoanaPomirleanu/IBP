# rescueTSS
# Decide which _0 TSS candidates to rescue based on promoter motifs.
# Only positionally conserved motifs are used, and a hit only counts at the
# position where MEME found the motif in the _1 promoters (motif_selection.tsv).
# Motif p-value cutoffs are calibrated on random genomic windows (scanBackground)
# with the same position rule, so the false-positive rate of the rescue as a
# whole (any motif) is measured, not assumed.
# Output: final table with a tier column (1 = _1 TSS, 2 = rescued _0, 3 = _0 not rescued).

RESCUE = config.get("motif rescue", {})
_RS_KM = config.get("known motifs") or {}
_RS_KM_ON = bool(_RS_KM.get("motif file") or _RS_KM.get("consensus"))


def _rescue_known_inputs(wc):
    """known-motif inputs, only when known motifs are given"""
    if not _RS_KM_ON:
        return {}
    d = f"results/promoter_motifs/{wc.sample}_{wc.ident}/motifs"
    return {"known_list": d + "/known_motifs.tsv",
            "training_fimo": d + "/fimo_training.tsv",
            "training_ids": d + "/training_ids.tsv"}


def _rescue_meme_file(wc):
    """file with the de novo training sites: meme.txt, or the combined file
    when only known motifs are used"""
    d = f"results/promoter_motifs/{wc.sample}_{wc.ident}"
    if _RS_KM_ON and not _RS_KM.get("use de novo motifs", True):
        return d + "/motifs/scan_motifs.meme"
    return d + "/meme/meme.txt"


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
    if _RS_KM_ON:
        d = f"results/promoter_motifs/{wc.sample}_{wc.ident}/motifs"
        o.append(f"--known-list {d}/known_motifs.tsv --training-fimo {d}/fimo_training.tsv "
                 f"--training-ids {d}/training_ids.tsv")
        for name, rng in (_RS_KM.get("positions") or {}).items():
            o.append(f"--known-position '{name}' {rng[0]} {rng[1]}")
    return " ".join(o)


rule rescueTSS:
    input:
        unpack(_rescue_known_inputs),                      # known motifs only
        random_fimo = "results/promoter_motifs/{sample}_{ident}/background/fimo_random.tsv",
        random_fa   = "results/promoter_motifs/{sample}_{ident}/background/random.fa",
        cand0_fimo  = rules.scanCandidates.output.fimo,  # takes whatever file scanCandidates writes as it's fimo output
        cand0_bed   = rules.selectCand0.output.bed,
        tss_fa      = SC_TSS_DIR + "/TSS_seq_{sample}_enriched_{ident}.fa.out",
        motifs      = _rescue_meme_file,                   # training sites + E-values
        scan_motifs = SC_MOTIF_FILE,                       # the file both scans used
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
        "logs/rescueTSS/{sample}_{ident}.log"
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
