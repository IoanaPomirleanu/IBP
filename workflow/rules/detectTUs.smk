# Transcription unit (TU) detection.
# A TU runs from a TSS to a downstream TTS on the same strand and is kept only if reads start at the TSS and end at the TTS.
# One TSS can give several TUs (read-through). A TSS whose reads end at no TTS gets one TU up to its furthest read 3' end.
# All positions below are 1-based; "points" are 1-nt BED intervals (pos-1, pos).

# 1. one point per TSS above the spike-in threshold (_1). TSS window: + [pos-up, pos+down], - [pos-down, pos+up]
rule tssPoints:
    input:
        'results/transcript_boundaries/TSS_{sample}/TSS_{sample}_{ident}/{sample}_enriched_{ident}_peaks.bed'
    output:
        'results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TSS_points_{sample}_{ident}.bed'
    params:
        down=config["TSS sequence extraction"]["downstream"]
    shell:
        """
        awk -F '\t' -v OFS='\t' -v d={params.down} '
            $4 ~ /_1$/ {{
                p = ($6 == "+") ? $3 - d : $2 + d
                print $1, p - 1, p, $4, 0, $6
            }}' {input} | sort -k1,1 -k2,2n > {output}
        """

# 2. one point per TTS. TTS window: + [pos-up, pos+down], - [pos-down, pos+up]; use the end, the start may be clamped to 0
rule ttsPoints:
    input:
        'results/transcript_boundaries/TTS_{sample}/TTS_{sample}_{ident}/TTS_{sample}_{ident}_corrected.bed'
    output:
        'results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TTS_points_{sample}_{ident}.bed'
    params:
        up=config["TTS sequence extraction"]["upstream"],
        down=config["TTS sequence extraction"]["downstream"]
    shell:
        """
        awk -F '\t' -v OFS='\t' -v u={params.up} -v d={params.down} '{{
                p = ($6 == "+") ? $3 - d : $3 - u
                print $1, p - 1, p, $4, 0, $6
            }}' {input} | sort -k1,1 -k2,2n > {output}
        """

# 3. every primary read of the enriched sample as a point at its 5' end, with its 3' end in the name column
rule readEnds:
    input:
        'results/alignments/BAM_files_{sample}/{sample}_enriched_{ident}.sorted.bam'
    output:
        temp('results/transcription_units/TU_{sample}/TU_{sample}_{ident}/read_ends_{sample}_{ident}.bed')
    conda:
        "../envs/env_tu.yaml"
    shell:
        """
        samtools view -b -F 0x904 {input} | bedtools bamtobed -i stdin | \
        awk -F '\t' -v OFS='\t' '{{
                if ($6 == "+") {{ f = $2 + 1; t = $3 }} else {{ f = $3; t = $2 + 1 }}
                print $1, f - 1, f, t, 0, $6
            }}' | sort -k1,1 -k2,2n > {output}
        """

# 4. reads starting within 'TSS window' nt of a TSS: a point at the read 3' end, name = TSS id, score = TSS position
rule readsPerTSS:
    input:
        tss='results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TSS_points_{sample}_{ident}.bed',
        reads='results/transcription_units/TU_{sample}/TU_{sample}_{ident}/read_ends_{sample}_{ident}.bed'
    output:
        'results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TSS_read_ends_{sample}_{ident}.bed'
    params:
        w=config["TU detection"]["TSS window"]
    conda:
        "../envs/env_tu.yaml"
    shell:
        """
        bedtools window -a {input.tss} -b {input.reads} -w {params.w} -sm | \
        awk -F '\t' -v OFS='\t' '{{ print $1, $10 - 1, $10, $4, $3, $6 }}' | \
        sort -k1,1 -k2,2n > {output}
        """

# 5. read 3' ends within 'TTS window' nt of a downstream TTS: count reads per (TSS, TTS), keep pairs with >= 'minimum reads'
# columns: chr, TSS id, TSS pos, TTS id, TTS pos, strand, n reads, type
rule pairTSSandTTS:
    input:
        ends='results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TSS_read_ends_{sample}_{ident}.bed',
        tts='results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TTS_points_{sample}_{ident}.bed'
    output:
        'results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TU_pairs_{sample}_{ident}.tsv'
    params:
        w=config["TU detection"]["TTS window"],
        m=config["TU detection"]["minimum reads"]
    conda:
        "../envs/env_tu.yaml"
    shell:
        """
        bedtools window -a {input.ends} -b {input.tts} -w {params.w} -sm | \
        awk -F '\t' -v OFS='\t' -v m={params.m} '
            ($6 == "+" && $9 > $5) || ($6 == "-" && $9 < $5) {{
                k = $4 SUBSEP $10
                n[k]++
                line[k] = $1 OFS $4 OFS $5 OFS $10 OFS $9 OFS $6
            }}
            END {{ for (k in n) if (n[k] >= m) print line[k], n[k], "TSS_TTS" }}' | \
        sort -k1,1 -k3,3n > {output}
        """

# 6. TSSs with reads but no TSS-TTS pair: TU up to the furthest read 3' end (max on +, min on -). Same columns as step 5
rule fallbackTUs:
    input:
        pairs='results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TU_pairs_{sample}_{ident}.tsv',
        ends='results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TSS_read_ends_{sample}_{ident}.bed'
    output:
        'results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TU_longest_read_{sample}_{ident}.tsv'
    shell:
        """
        awk -F '\t' -v OFS='\t' '
            FILENAME == ARGV[1] {{ paired[$2] = 1; next }}
            !($4 in paired) {{
                k = $4; n[k]++
                if (!(k in e) || ($6 == "+" && $3 > e[k]) || ($6 == "-" && $3 < e[k])) e[k] = $3
                c[k] = $1; t[k] = $5; s[k] = $6
            }}
            END {{ for (k in n) print c[k], k, t[k], "longest_read", e[k], s[k], n[k], "longest_read" }}
            ' {input.pairs} {input.ends} | sort -k1,1 -k3,3n > {output}
        """

# 7. all TUs as BED6 (score = supporting reads, for IGV) and as a table with header
rule combineTUs:
    input:
        'results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TU_pairs_{sample}_{ident}.tsv',
        'results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TU_longest_read_{sample}_{ident}.tsv'
    output:
        bed='results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TU_{sample}_{ident}.bed',
        tsv='results/transcription_units/TU_{sample}/TU_{sample}_{ident}/TU_{sample}_{ident}.tsv'
    shell:
        """
        cat {input} | sort -k1,1 -k3,3n | awk -F '\t' -v OFS='\t' -v tsv={output.tsv} '
            BEGIN {{ print "TU_id", "chr", "start", "end", "strand", "TSS", "TSS_pos", "end_id", "end_pos", "n_reads", "type" > tsv }}
            {{
                a = ($3 < $5) ? $3 : $5
                b = ($3 < $5) ? $5 : $3
                id = "TU_" NR
                print $1, a - 1, b, id, $7, $6
                print id, $1, a - 1, b, $6, $2, $3, $4, $5, $7, $8 > tsv
            }}' > {output.bed}
        """
