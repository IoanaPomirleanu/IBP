# =============================================================================
#  Promoter motif discovery (per {sample}_{ident})
#
#  dedup_training         : collapse redundant _1 training windows
#  motif_background_model : composition model of the phage genome
#  discover_motifs        : MEME on the deduplicated _1 windows
#  (candidate selection and FIMO scanning: rules/scan_candidates.smk)
#
#  Input : results/promoter_motifs/{sample}_{ident}/meme_training.fa
#          (_1 windows, rule selectAboveThresholdSequences)
#  Output: results/promoter_motifs/{sample}_{ident}/meme/meme.txt
#          = the default "motif file" of scan_candidates (and scan_background)
#
#  Include from workflow/Snakefile:   include: "rules/motifs.smk"
#  Target: motifs_all (all samples in the config), or the meme.txt path
# =============================================================================
import os

MOT = config.get("motifs", {})
MOT_DIR = "results/promoter_motifs/{sample}_{ident}"
MOT_SCRIPTS = os.path.join(workflow.basedir, "scripts")
MOT_ORDER = MOT.get("markov order", 2)


rule motifs_all:
    input:
        expand(MOT_DIR + "/meme/meme.txt",
               sample=config["Sample name"], ident=config["ID"]),


# ---------------------------------------------------------------------------
# training windows
# ---------------------------------------------------------------------------
rule dedup_training:
    input:
        MOT_DIR + "/meme_training.fa",
    output:
        MOT_DIR + "/discovery/meme_training_dedup.fa",
    params:
        script=os.path.join(MOT_SCRIPTS, "dedup_windows.py"),
        max_offset=MOT.get("dedup distance") or config["cluster width"]["TSS"],
    log:
        "logs/motifs/{sample}_{ident}/dedup_windows.log",
    conda:
        "../envs/env_motif.yaml"
    shell:
        "python {params.script} --fasta {input} --out {output} "
        "--max-offset {params.max_offset} 2> {log}"


# Background model of the phage genome (base and short-word composition).
# MEME uses it to judge how surprising a pattern is in THIS genome.
# Genome: "genome fasta" (default: "fasta file"). Records: "genome records"
# (default: the records named in the training headers, e.g. LUZ7), so a
# concatenated host + phage + spike-in FASTA is reduced to the phage.
rule motif_background_model:
    input:
        genome=MOT.get("genome fasta") or config["fasta file"],
        training=MOT_DIR + "/meme_training.fa",
    output:
        genome=MOT_DIR + "/discovery/background_genome.fa",
        model=MOT_DIR + f"/discovery/background_order{MOT_ORDER}.txt",
    params:
        script=os.path.join(MOT_SCRIPTS, "extract_records.py"),
        order=MOT_ORDER,
        records=("--records " + " ".join(MOT["genome records"])) if MOT.get("genome records") else "",
    log:
        "logs/motifs/{sample}_{ident}/background_model.log",
    conda:
        "../envs/env_motif.yaml"
    shell:
        """
        python {params.script} {input.genome} {output.genome} {params.records} \
            --records-from {input.training} > {log} 2>&1
        fasta-get-markov -m {params.order} -dna {output.genome} {output.model} >> {log} 2>&1
        """


rule discover_motifs:
    input:
        training=rules.dedup_training.output[0],
        bfile=rules.motif_background_model.output.model,
    output:
        MOT_DIR + "/meme/meme.txt",
    params:
        outdir=MOT_DIR + "/meme",
        min_seqs=MOT.get("min sequences", 5),
        nmotifs=MOT.get("nmotifs", 3),
        minw=MOT.get("minw", 6),
        maxw=MOT.get("maxw", 20),
        mode=MOT.get("mode", "zoops"),
    log:
        "logs/motifs/{sample}_{ident}/meme.log",
    conda:
        "../envs/env_motif.yaml"
    shell:
        """
        n=$(grep -c '>' {input.training} || true)
        if [ "$n" -ge {params.min_seqs} ]; then
            # windows are already oriented 5'->3' along the transcript,
            # so MEME searches the given strand only (no -revcomp)
            meme {input.training} -dna -oc {params.outdir} -mod {params.mode} \
                 -nmotifs {params.nmotifs} -minw {params.minw} -maxw {params.maxw} \
                 -bfile {input.bfile} -nostatus > {log} 2>&1
        else
            mkdir -p {params.outdir}
            printf 'MEME version 5\\n' > {output}
            echo "only $n independent windows (< {params.min_seqs}); discovery skipped" > {log}
        fi
        """
