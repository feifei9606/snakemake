rule star_index:
    input:
        fasta = config["ref"]["genome"],
        ann   = config["ref"]["annotation"]
    output:
        done = touch("01_data/index/star/makeidx.done")
    threads: 16
    log: "logs/index/star_index.log"
    conda : "../envs/circ.yaml"
    shell : 
        """
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir "01_data/index/star" --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.ann} 2> {log}
        """

rule star_align:
    input:
        idxdone = "01_data/index/star/makeidx.done",
        sample = get_trimmed_reads
    output:
        "04_align/star/{sample}/Aligned.sortedByCoord.out.bam",
        "04_align/star/{sample}/Chimeric.out.junction"
    log: "logs/align/star/{sample}.log"
    params:
        index="01_data/index/star",
        ann = config["ref"]["annotation"],
        prefix = "04_align/star/{sample}/"
    threads: 16
    conda : "../envs/circ.yaml"
    shell :
        """
            STAR --runThreadN {threads} --genomeDir {params.index} --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate \
            --readFilesIn {input.sample} --readFilesCommand zcat --outReadsUnmapped Fastx --outSJfilterOverhangMin 15 15 15 15 \
            --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --outFilterMultimapNmax 20 --outFilterScoreMin 1 --outFilterMatchNmin 1 \
            --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --sjdbGTFfile {params.ann} --outStd Log >{log}
        """

rule star_align_mate:
    input:
        idxdone = "01_data/index/star/makeidx.done",
        sample = "03_trim/{sample}.{mate}.trim.fastq.gz"
    output:
        "04_align/star/{sample}/{mate}/Chimeric.out.junction"
    log: "logs/align/star/{sample}.{mate}.log"
    params:
        index="01_data/index/star",
        ann = config["ref"]["annotation"],
        prefix = "04_align/star/{sample}/{mate}/"
    threads: 8
    conda : "../envs/circ.yaml"
    shell :
        """
            STAR --runThreadN {threads} \
            --genomeDir {params.index} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype None \
            --readFilesIn {input.sample} \
            --readFilesCommand zcat \
            --outReadsUnmapped Fastx \
            --outSJfilterOverhangMin 15 15 15 15 \
            --alignSJoverhangMin 15 \
           --alignSJDBoverhangMin 15 \
           --seedSearchStartLmax 30 \
           --outFilterMultimapNmax 20 \
           --outFilterScoreMin 1 \
           --outFilterMatchNmin 1 \
           --outFilterMismatchNmax 2 \
           --chimSegmentMin 15 \
           --chimScoreMin 15 \
           --chimScoreSeparation 10 \
           --chimJunctionOverhangMin 15 \
           --outStd Log >{log}
        """

rule bwa_index:
    input:
        config["ref"]["genome"]
    output:
        touch("01_data/index/bwa/makeidx.done")
    log:
        "logs/index/bwa_index.log"
    params:
        prefix="01_data/index/bwa/" + config["ref"]["prefix"],
        algorithm="bwtsw"
    conda : "../envs/circ.yaml"
    shell : 
        """
        bwa index -a {params.algorithm} -p {params.prefix} {input} 2> {log}
        """

rule bwa_align:
    input: 
        samples = get_trimmed_reads,
        idxdone = "01_data/index/bwa/makeidx.done"
    output: "04_align/bwa/{sample}_bwa.sam"
    log: "logs/align/bwa/{sample}.log"
    params:
        index="01_data/index/bwa/" + config["ref"]["prefix"]
    threads: 16
    conda : "../envs/circ.yaml"
    shell:
        """
        bwa mem -t {threads} -T 19 {params.index} {input.samples} > {output} 2> {log}
        """