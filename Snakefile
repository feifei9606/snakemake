###############################################################################################################################################
## Configuration file
###############################################################################################################################################
import os,sys
import pandas as pd
if len(config) == 0:
    if os.path.isfile("config.yaml"):
        configfile : "config.yaml"
    else:
        sys.exit("Make sure there is a config.yaml in " + os.getcwd() + " or specify one with the --configfile commandline parameter.")


genome = config["refs"]["genome"]
transcriptome_gtf = config["refs"]["transcriptome_gtf"]
indexpath = config["index"]["indexpath"]
indexname = config["index"]["indexname"]
###############################################################################################################################################
## Sample and conditions
###############################################################################################################################################

units = pd.read_csv(config["units"], dtype=str,sep = "\t").set_index(["sample"], drop=False)
samplefile = config["units"]

SAMPLES = units.index.get_level_values("sample").unique().tolist()

###############################################################################################################################################
## Input functions for rules
###############################################################################################################################################
def get_sra(wildcards):
    return units.loc[(wildcards.sample),["sra"]].dropna()

def get_gtf(wildcards):
    return expand("mapped/gtf/{sample}.gtf",sample=SAMPLES)
###############################################################################################################################################
## Desired outputs expand("mapped/gtf/{sample}.gtf",sample=SAMPLES)
###############################################################################################################################################
rule all:
    input:
        # FASTQC_BEFORE_TRIMMING = expand("fastqc/fastqc_before_trimming/{sample}_fastqc.html", sample = SAMPLES),
        # FASTQC_AFTER_TRIMMING = expand("fastqc/fastqc_before_trimming/{sample}_fastqc.html", sample = SAMPLES),
        # GTF = "mapped/gtf_merged/merged.gtf",
        # COUNTS = "mapped/count/readcount/count.txt",
        DESeq2 = "Rplot/DESeq2_results_all.csv"
    message:
        "Job Done! Removing temportory directory."
###############################################################################################################################################
## sra2fastq
###############################################################################################################################################

rule sra2fastq:
    input:
        sra = get_sra,
    output:
        fastq = "fastq/{sample}.fastq"
    message:
        "{wildcards.sample}.sra to {wildcards.sample}.fastq"
    log:
        "logs/sra2fastq/{sample}.log"
    params:
        outdir = "fastq"
    shell:
        "fastq-dump -O {params.outdir} {input} 1>{log} 2>&1"

###############################################################################################################################################
## Fastq QC before / after trimming
###############################################################################################################################################

rule fastqc_before_trimming:
    input:
        fq = "fastq/{sample}.fastq"
    output:
        html = "fastqc/fastqc_before_trimming/{sample}_fastqc.html",
        zipfile = "fastqc/fastqc_before_trimming/{sample}_fastqc.zip"
    message:
        "Quality check of {wildcards.sample} sample with fastqc"
    log:
        "logs/fastqc_before_trimming/{sample}.log"
    params:
        outdir = "fastqc/fastqc_before_trimming",
        threads = 4
    shell:
        "fastqc -f fastq -o {params.outdir} -t {params.threads} {input} 1>{log} 2>&1"

rule qc_trimmer:
    input:
        fq = "fastq/{sample}.fastq"
    output:
        trimmedfq = "trimmed/fastx_trimmer/{sample}_trimmed.fastq"
    message:
        "Trimming {wildcards.sample} reads"
    log:
        "logs/trimmed_fastx_trimmer/{sample}.log"
    shell:
        "fastx_trimmer -Q 33 -f 12 -i {input} -o {output} 1>{log} 2>&1"

rule qc_filter:
    input:
        trimmedfq = "trimmed/fastx_trimmer/{sample}_trimmed.fastq"
    output:
        filterdfq = "trimmed/fastq_quality_filter/{sample}_filtered.fastq"
    message:
        "Filtering {wildcards.sample} reads"
    log:
        "logs/trimmed_fastq_quality_filter/{sample}.log"
    shell:
        "fastq_quality_filter -Q 33 -q 20 -p 80 -i {input} -o {output} 1>{log} 2>&1"

rule fastqc_after_trimming:
    input:
        fastq = "trimmed/fastq_quality_filter/{sample}_filtered.fastq"
    output:
        html = "fastqc/fastqc_after_qc/{sample}_filtered_fastqc.html",
        zipfile = "fastqc/fastqc_after_qc/{sample}_filtered_fastqc.zip"
    message:
        "Quality check of filterd {wildcards.sample} sample with fastqc"
    log:
        "logs/fastqc_after_qc/{sample}.log"
    params:
        outdir = "fastqc/fastqc_after_qc",
        threads = 4
    shell:
        "fastqc -o {params.outdir} -t {params.threads} {input} 1>{log} 2>&1"

###############################################################################################################################################
## RNA-seq read alignment
###############################################################################################################################################

rule buildindex:
    input:
        genome = genome,
        transcriptome_gtf = transcriptome_gtf
    output:
        [ "{}/{}.{}.ht2".format(indexpath,indexname,i) for i in range(1,9)]
    message:
        "Build HiSAT2 genome index"
    params:
        indexpath = indexpath,
        indexname = indexname
    log:
        "logs/buildindex/hisat2_buildindex.log"
    threads : 20
    shell:
        "hisat2_extract_splice_sites.py {input.transcriptome_gtf} > {input.transcriptome_gtf}.ss |"
        "hisat2_extract_exons.py {input.transcriptome_gtf} > {input.transcriptome_gtf}.exon |"
        "hisat2-build -p {threads} --ss {input.transcriptome_gtf}.ss --exon {input.transcriptome_gtf}.exon {input.genome} {params.indexpath}/{params.indexname} 1>{log} 2>&1"


rule hisat_mapping:
    input:
        fq = "trimmed/fastq_quality_filter/{sample}_filtered.fastq",
        indexfiles = [ "{}/{}.{}.ht2".format(indexpath,indexname,i) for i in range(1,9)]
    output:
        sam = "mapped/sam/{sample}.sam"
    message:
        "Mapping reads to genome to sam files."
    log:
        "logs/mapped/{sample}.log"
    params:
        index = indexpath + "/" + indexname
    threads: 10
    shell:
        "hisat2 -p {threads} --dta -x {params.index} -U {input.fq} -S {output} 1>{log} 2>&1"

rule sam2bam:
    input:
        sam = "mapped/sam/{sample}.sam"
    output:
        bam = "mapped/bam/{sample}.bam"
    message:
        "Convert sam To bam"
    log:
        "logs/sam2bam/{sample}.log"
    threads: 5
    shell:
        "samtools sort -@ {threads} -m 1000M -o {output.bam} {input} 1>{log} 2>&1"

rule stringtie:
    input:
        transcriptome_gtf = transcriptome_gtf,
        bam = "mapped/bam/{sample}.bam"
    output:
        sample_gtf = "mapped/gtf/{sample}.gtf"
    message:
        "creating transcriptome to stringtie_transcriptome.gtf"
    log:
        "logs/stringtie/{sample}.log"
    threads: 10
    shell:
        "stringtie -p {threads} -G {input.transcriptome_gtf} -o {output.sample_gtf} -l {wildcards.sample} {input.bam} 1>{log} 2>&1 |"
        "echo {output.sample_gtf} >> mapped/gtf/gtflist.txt "

rule stringtie_merge:
    input:
        sample_gtf = get_gtf,
        transcriptome_gtf = transcriptome_gtf
    output:
        merged_gtf = "mapped/gtf_merged/merged.gtf"
    message:
        "Merge gtf generated by stringtie"
    log:
        "logs/stringtie/merged.log"
    threads: 10
    shell:
        "stringtie --merge -p {threads} -G {input.transcriptome_gtf}  -o {output.merged_gtf} mapped/gtf/gtflist.txt 1>{log} 2>&1"

rule calfpkm:
    input:
        merged_gtf = "mapped/gtf_merged/merged.gtf",
        bam = "mapped/bam/{sample}.bam"
    output:
        sample_genes_fpkm = "mapped/count/fpkm/{sample}_genes.gtf",
        sample_transcripts_fpkm = "mapped/count/fpkm/{sample}_transcripts.gtf"
    message:
        "Calculate fpkm using stringtie"
    log:
        "logs/calfpkm/{sample}.log"
    threads: 10
    shell:
        "stringtie -e -p {threads} -G {input.merged_gtf} -A {output.sample_genes_fpkm} -o {output.sample_transcripts_fpkm} {input.bam} 1>{log} 2>&1"

rule calreadcount:
    input:
        merged_gtf = "mapped/gtf_merged/merged.gtf",
        bams = expand("mapped/bam/{sample}.bam",sample=SAMPLES)
    output:
        total_count = "mapped/count/readcount/count.txt",
    message:
        "Calculate reads count using featureCounts"
    log:
        "logs/calreadcount/totalcount.log"
    threads: 10
    shell:
        "featureCounts -T {threads} -O -t transcript -g transcript_id -a {input.merged_gtf} -o {output.total_count} {input.bams} 1>{log} 2>&1"


###############################################################################################################################################
## RNA-seq downstream analysis
###############################################################################################################################################

rule DESeq2_analysis:
    input:
        counts = "mapped/count/readcount/count.txt",
        samplefile = samplefile,
    output:
        DESeq2_results_all = "Rplot/DESeq2_results_all.csv"
    message:
        "Calculate the differential expression using DESeq2"
    log:
        "logs/DESeq2/DESeq2.log"
    shell:
        "Rscript scripts/DESeq2.R -c {input.counts} -s {input.samplefile} -o {output.DESeq2_results_all} 1>{log} 2>&1 "