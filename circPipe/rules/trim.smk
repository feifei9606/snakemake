wildcard_constraints:
    sample="[a-zA-Z0-9]+"

rule trim_fastp_se:
    input : sample=["01_data/raw/{sample}.fastq.gz"]
    output : 
        trimmed="03_trim/{sample}.trim.fastq.gz",
        html="logs/trim/{sample}.html",
        json="logs/trim/{sample}.json"
    threads : 8
    log:
        "logs/trim/{sample}_se.log"
    params:
        extra=""
    wrapper: "0.45.1/bio/fastp"

rule trim_fastp_pe:
    input:
        sample=["01_data/raw/{sample}.R1.fastq.gz","01_data/raw/{sample}.R2.fastq.gz"]
    output:
        trimmed = ["03_trim/{sample}.R1.trim.fastq.gz","03_trim/{sample}.R2.trim.fastq.gz"],
        html = "logs/trim/{sample}.html",
        json = "logs/trim/{sample}.json"
    threads: 8
    params:
        extra=""
    log:
        "logs/trim/{sample}_pe.log"
    wrapper:
        "0.45.1/bio/fastp"