"""
scanBackground
Null model for calibrating the motif threshold
FIMO is run on random genomic windows
-> how often motifs hit in real, non-promoter DNA
rescueTSS uses these hits to choose a p-value with a known false-positive rate
"""

#Zero-order background model (base frequencies of the genome)
# scanCandidates should use this same file (--bgfile)
rule backgroundModel: 
    input: 
        genome= rules.motifBackgroundModel.output.genome   # phage-only genome (motifs.smk)

    output:
        bg = "results/promoter_motifs/{sample}_{ident}/bg.txt"
    log:
        "logs/backgroundModel/{sample}_{ident}.log"
    conda:
        "../envs/env_scan_background.yaml"
    shell:
        "fasta-get-markov -m 0 {input.genome} {output.bg} 2> {log}"


rule scanBackground: 
    input: 
        genome   = rules.motifBackgroundModel.output.genome,   # phage-only genome (motifs.smk)
        windows1 = "results/promoter_motifs/{sample}_{ident}/meme_training.fa",  # only used for the window length
        motifs   =  SC_MOTIF_FILE,
        bg       = "results/promoter_motifs/{sample}_{ident}/bg.txt"

    output:
        sizes       = "results/promoter_motifs/{sample}_{ident}/background/genome.sizes",
        random_bed  = "results/promoter_motifs/{sample}_{ident}/background/random.bed",
        random_fa   = "results/promoter_motifs/{sample}_{ident}/background/random.fa",
        fimo_random = "results/promoter_motifs/{sample}_{ident}/background/fimo_random.tsv"

    params:
        n_random  = config["motif background"]["n random windows"],
        thresh    = config["motif background"]["fimo threshold"],
        fimo_opts = config["motif background"]["fimo options"],
        seed      = 42
    
    log: 
        "logs/scanBackground/{sample}_{ident}.log"
    conda:
        "../envs/env_scan_background.yaml"
    shell: 
         """
        (
        # window length = length of peaks_1 windows, so the random windows match
        W=$(awk '!/^>/ {{print length($0); exit}}' {input.windows1})
        echo "window length: $W"

        # chromosome names + lengths, needed by bedtools random
        samtools faidx {input.genome}
        cut -f1,2 {input.genome}.fai > {output.sizes}

        # random genomic windows on random strands, sequence taken strand-aware (-s)
        # random windows are made in two steps: choose where they are => fetch their seq
        bedtools random -l $W -n {params.n_random} -seed {params.seed} \
            -g {output.sizes} > {output.random_bed}
        bedtools getfasta -fi {input.genome} -bed {output.random_bed} -s -nameOnly \
            -fo {output.random_fa}

        # FIMO on the random windows (same settings as scanCandidates)
        fimo --text {params.fimo_opts} --thresh {params.thresh} --bgfile {input.bg} \
            {input.motifs} {output.random_fa} > {output.fimo_random}
        ) 2> {log}
        """