# promoter windows of candidates at or above the spike-in threshold (_1), on the phage contig only
rule selectAboveThresholdSequences:
    input: 'results/transcript_boundaries/TSS_{sample}/TSS_{sample}_{ident}/TSS_seq_{sample}_enriched_{ident}.fa.out'
    output: 'results/promoter_motifs/{sample}_{ident}/meme_training.fa'
    shell:
        """
        awk -v phage="{wildcards.sample}" '
          /^>/ {{ split(substr($0, 2), a, "::"); keep = (a[1] ~ /_1$/ && index(a[2], phage) > 0) }}
          keep
        ' {input} > {output}
        """
