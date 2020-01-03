rule before_fastqc_se:
    input:
        "01_data/raw/{sample}.fastq.gz"
    output:
        html="02_qc/before_trimmed/{sample}_fastqc.html",
        zip="02_qc/before_trimmed/{sample}_fastqc.zip"
    log:
        "logs/before_fastqc/{sample}_se.log"
    conda : "../envs/circ.yaml"
    shell : 
        """
        fastqc {input} -o 02_qc/before_trimmed 2> {log}
        """
rule before_fastqc_pe:
    input:
        r1 = "01_data/raw/{sample}.R1.fastq.gz",
        r2 = "01_data/raw/{sample}.R2.fastq.gz"
    output:
        html1="02_qc/before_trimmed/{sample}.R1_fastqc.html",
        zip1="02_qc/before_trimmed/{sample}.R1_fastqc.zip",
        html2="02_qc/before_trimmed/{sample}.R2_fastqc.html",
        zip2="02_qc/before_trimmed/{sample}.R2_fastqc.zip"
    log:
        "logs/before_fastqc/{sample}_pe.log"
    threads: 2
    conda : "../envs/circ.yaml"
    shell : 
        """
        fastqc {input.r1} {input.r2} -o 02_qc/before_trimmed -t {threads} 2> {log}
        """
rule after_fastqc:
    input:
        "03_trim/{sample}.trim.fastq.gz"
    output:
        html="02_qc/after_trimmed/{sample}.trim_fastqc.html",
        zip="02_qc/after_trimmed/{sample}.trim_fastqc.zip"
    log:
        "logs/after_fastqc/{sample}_se.log"
    conda : "../envs/circ.yaml"
    shell : 
        """
        fastqc {input} -o 02_qc/after_trimmed 2> {log}
        """
rule after_fastqc_pe:
    input:
        r1 = "03_trim/{sample}.R1.trim.fastq.gz",
        r2 = "03_trim/{sample}.R2.trim.fastq.gz"
    output:
        html1="02_qc/after_trimmed/{sample}.R1.trim_fastqc.html",
        zip1="02_qc/after_trimmed/{sample}.R1.trim_fastqc.zip",
        html2="02_qc/after_trimmed/{sample}.R2.trim_fastqc.html",
        zip2="02_qc/after_trimmed/{sample}.R2.trim_fastqc.zip"
    log:
        "logs/after_fastqc/{sample}_pe.log"
    threads : 2
    conda : "../envs/circ.yaml"
    shell : 
        """
        fastqc {input.r1} {input.r2} -o 02_qc/after_trimmed -t {threads} 2> {log}
        """
rule multiqc_se:
    input:
        expand(["02_qc/before_trimmed/{u.sample}_fastqc.zip","02_qc/after_trimmed/{u.sample}.trim_fastqc.zip"],u=units.itertuples())
    output:
        "02_qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.45.1/bio/multiqc"

rule multiqc_pe:
    input: 
        expand(["02_qc/before_trimmed/{u.sample}.{read}_fastqc.zip"],read=("R1","R2"),u=units.itertuples()),
        expand(["02_qc/after_trimmed/{u.sample}.{read}.trim_fastqc.zip"],read=("R1","R2"),u=units.itertuples()),
        expand(["04_align/star/{u.sample}"],u=units.itertuples())
    output:
        "02_qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.45.1/bio/multiqc"

