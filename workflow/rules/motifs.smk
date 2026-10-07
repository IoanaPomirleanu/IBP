# =============================================================================
#  Promoter motif module
#
#  dedup_training         : collapse redundant _1 training windows
#  motif_background_model : composition model of the phage genome
#  discover_motifs        : MEME on the deduplicated _1 windows
#  (candidate selection and FIMO scanning: rules/scan_candidates.smk)
#
#  Include from workflow/Snakefile:   include: "rules/motifs.smk"
#  Targets:  discover_motifs / motifs_all   (motif discovery)
# =============================================================================
import os

MOT = config["motifs"]
SAMPLE = config["Sample name"]
MOT_OUT = f"{MOT.get('outdir', 'results/motifs')}/{SAMPLE}"
MOT_SCRIPTS = os.path.join(workflow.basedir, "scripts")

rule motifs_all:
    input:
        f"{MOT_OUT}/meme/meme.txt",


# ---------------------------------------------------------------------------
# training windows
# ---------------------------------------------------------------------------
rule dedup_training:
    input:
        MOT["training fasta"],                 # meme_training.fa (_1 windows)
    output:
        f"{MOT_OUT}/meme_training_dedup.fa",
    params:
        script=os.path.join(MOT_SCRIPTS, "dedup_windows.py"),
        max_offset=MOT.get("dedup distance") or config["cluster width"]["TSS"],
    log:
        f"{MOT_OUT}/logs/dedup_windows.log",
    conda:
        "../envs/env_motif.yaml"
    shell:
        "python {params.script} --fasta {input} --out {output} "
        "--max-offset {params.max_offset} 2> {log}"


def _background_records():
    """Genome records for the background model.

    "genome records" in the config wins; otherwise the records named in the
    training window headers (name::RECORD:start-end(strand)).
    """
    if MOT.get("genome records"):
        return list(MOT["genome records"])
    recs = set()
    try:
        with open(MOT["training fasta"]) as fh:
            for line in fh:
                if line.startswith(">") and "::" in line:
                    recs.add(line[1:].split("::", 1)[1].rsplit(":", 1)[0])
    except FileNotFoundError:
        pass
    return sorted(recs)


MOT_RECORDS = _background_records()


# Background model of the phage genome (base and short-word composition).
# MEME and FIMO use it to judge how surprising a match is in THIS genome.
rule motif_background_model:
    input:
        MOT.get("genome fasta") or config["fasta file"],
    output:
        genome=f"{MOT_OUT}/background_genome.fa",
        model=f"{MOT_OUT}/background_order{MOT.get('markov order', 2)}.txt",
    params:
        script=os.path.join(MOT_SCRIPTS, "extract_records.py"),
        order=MOT.get("markov order", 2),
        records=" ".join(MOT_RECORDS),
    log:
        f"{MOT_OUT}/logs/background_model.log",
    conda:
        "../envs/env_motif.yaml"
    shell:
        """
        python {params.script} {input} {output.genome} {params.records} > {log} 2>&1
        fasta-get-markov -m {params.order} -dna {output.genome} {output.model} >> {log} 2>&1
        """


rule discover_motifs:
    input:
        training=rules.dedup_training.output[0],
        bfile=rules.motif_background_model.output.model,
    output:
        f"{MOT_OUT}/meme/meme.txt",
    params:
        outdir=f"{MOT_OUT}/meme",
        min_seqs=MOT.get("min sequences", 5),
        nmotifs=MOT.get("nmotifs", 3),
        minw=MOT.get("minw", 6),
        maxw=MOT.get("maxw", 20),
        mode=MOT.get("mode", "zoops"),
    log:
        f"{MOT_OUT}/logs/meme.log",
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
