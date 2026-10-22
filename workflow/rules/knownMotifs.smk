# =============================================================================
#  Known (user-supplied) promoter motifs
#
#  combine_motifs : one motif file for both FIMO scans = de novo MEME motifs
#                   (optional) + known motifs (KNOWN-<name>), from a MEME-format
#                   file and/or IUPAC consensus sequences in the config
#  scan_training  : FIMO on the deduplicated _1 windows, so rescue_tss can see
#                   where each known motif sits relative to the TSS (and how many
#                   _1 promoters it hits)
#
#  Active only when "known motifs" has a "motif file" or "consensus" entry.
#  The motif file then used by scan_background and scan_candidates is
#  results/promoter_motifs/{sample}_{ident}/motifs/scan_motifs.meme, so the
#  known motifs get the same random-window calibration as the de novo ones.
# =============================================================================

KM = config.get("known motifs") or {}
KM_ON = bool(KM.get("motif file") or KM.get("consensus"))
KM_DENOVO = KM.get("use de novo motifs", True)
KM_DIR = "results/promoter_motifs/{sample}_{ident}/motifs"


def _km_inputs(wc):
    d = {}
    if KM_DENOVO:
        d["denovo"] = f"results/promoter_motifs/{wc.sample}_{wc.ident}/meme/meme.txt"
    if KM.get("motif file"):
        d["known"] = KM["motif file"]
    return d


def _km_args(wc, input):
    o = []
    if "denovo" in input.keys():
        o.append(f"--denovo {input.denovo}")
    if "known" in input.keys():
        o.append(f"--known-file {input.known}")
    cons = KM.get("consensus") or {}
    if cons:
        o.append("--consensus " + " ".join(f"'{k}={v}'" for k, v in cons.items()))
    o.append(f"--mismatch {KM.get('consensus mismatch', 0.05)}")
    return " ".join(o)


rule combine_motifs:
    input:
        unpack(_km_inputs),
    output:
        motifs=KM_DIR + "/scan_motifs.meme",
        known=KM_DIR + "/known_motifs.tsv",
    params:
        script=workflow.basedir + "/scripts/prepare_known_motifs.py",
        args=_km_args,
    log:
        "logs/combine_motifs/{sample}_{ident}.log",
    conda:
        "../envs/env_scan_background.yaml"
    shell:
        "python {params.script} {params.args} "
        "--out {output.motifs} --out-list {output.known} 2> {log}"


rule scan_training:
    input:
        fa="results/promoter_motifs/{sample}_{ident}/discovery/meme_training_dedup.fa",
        motifs=KM_DIR + "/scan_motifs.meme",
        bg="results/promoter_motifs/{sample}_{ident}/bg.txt",
    output:
        fimo=KM_DIR + "/fimo_training.tsv",
        ids=KM_DIR + "/training_ids.tsv",
    params:
        fa=KM_DIR + "/training_scan.fa",
        thresh=config["motif background"]["fimo threshold"],
        fimo_opts=config["motif background"]["fimo options"],
    log:
        "logs/scan_training/{sample}_{ident}.log",
    conda:
        "../envs/env_scan_background.yaml"
    shell:
        """
        (
        # short IDs (trn0001, ...) so FIMO cannot read the names as coordinates
        awk '/^>/ {{n++; id = sprintf("trn%04d", n); print id "\\t" substr($1, 2) > "{output.ids}"; print ">" id; next}} {{print}}' \
            {input.fa} > {params.fa}
        if grep -q '^MOTIF' {input.motifs} && [ -s {params.fa} ]; then
            fimo --text {params.fimo_opts} --thresh {params.thresh} --bgfile {input.bg} \
                {input.motifs} {params.fa} > {output.fimo}
        else
            printf 'motif_id\\tmotif_alt_id\\tsequence_name\\tstart\\tstop\\tstrand\\tscore\\tp-value\\tq-value\\tmatched_sequence\\n' > {output.fimo}
        fi
        touch {output.ids}
        ) 2> {log}
        """
