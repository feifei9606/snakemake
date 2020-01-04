junc_pattern = "04_align/star/{u.sample}/{mate}Chimeric.out.junction"

def get_sheets(wildcards):
    if not peFlag:
        return ["05_output/dcc/samplesheet_se"]
    return ["05_output/dcc/samplesheet_pe",
            "05_output/dcc/mate1",
            "05_output/dcc/mate2"]

def get_aligns(wc):
    return expand(
        "04_align/star/{u.sample}/Aligned.sortedByCoord.out.bam",
        u=units.itertuples())
def mates(wildcards,input):
    if peFlag:
        # paired-end sample
        return "-mt1 @{} -mt2 @{}".format(input.sheet[1],input.sheet[2])
    # single end sample
    return ""
def get_pi(wildcards):
    if peFlag:
        return "-Pi"
    return ""


rule ciri:
    input: "04_align/bwa/{sample}_bwa.sam"
    output: "05_output/ciri/{sample}_circRNA.csv"
    log : "logs/ciri/{sample}.log"
    params:
        ref = config["ref"]["genome"],
        ann = config["ref"]["annotation"]
    threads: 16
    shell : 
        """
        CIRI2.pl -I {input} -O {output} -F {params.ref} -A {params.ann} -T {threads} -G {log}
        """
rule dcc_se_samplesheet:
    input: expand(junc_pattern,u=units.itertuples(),mate="")
    output: "05_output/dcc/samplesheet_se"
    run :
        with open(output[0],"w") as f:
            f.write("\n".join(input))

rule dcc_pe_samplesheet:
    input: 
        sample = expand(junc_pattern,u=units.itertuples(),mate=""),
        mate1 = expand(junc_pattern,u=units.itertuples(),mate="R1/"),
        mate2 = expand(junc_pattern,u=units.itertuples(),mate="R2/")
    output:
        "05_output/dcc/samplesheet_pe",
        "05_output/dcc/mate1",
        "05_output/dcc/mate2"
    run :
        with open(output[0], "w") as f:
            f.write("\n".join(input.sample))
        with open(output[1], "w") as f:
            f.write("\n".join(input.mate1))
        with open(output[2], "w") as f:
            f.write("\n".join(input.mate2))

rule dcc:
    input:
        sheet = get_sheets,
        align = get_aligns
    output : "05_output/dcc/CircRNACount"
    params:
        ref = config["ref"]['genome'],
        ann = config["ref"]['annotation'],
        wkdir = "05_output/dcc",
        mates = mates,
        pi = get_pi
    threads : 16
    log : "logs/dcc/dcc.log"
    shell:
        """
        source activate py2;
        DCC @{input.sheet[0]} {params.mates} -F -D -fg -G -Nr 5 1 {params.pi} -T {threads} -A {params.ref} -an {params.ann} -O {params.wkdir} -t {params.wkdir}/_tmp 2>{log}
        """

rule circexplorer_parse_star:
    input: "04_align/star/{sample}/Chimeric.out.junction"
    output: "05_output/circexplorer/star/{sample}/bsj.bed"
    log:  "logs/circexplorer/star/{sample}_parse.log"
    conda : "../envs/circ.yaml"
    shell :
        """
            CIRCexplorer2 parse -t STAR -b {output} {input} > {log}
        """

rule circexplorer_parse_bwa :
    input: "04_align/bwa/{sample}_bwa.sam"
    output: "05_output/circexplorer/bwa/{sample}/bsj.bed"
    log:  "logs/circexplorer/bwa/{sample}_parse.log"
    conda : "../envs/circ.yaml"
    shell :
        """
            CIRCexplorer2 parse -t BWA -b {output} {input} > {log}
        """

rule circexplorer_annotate :
    input: 
        bsj = "05_output/circexplorer/{aligner}/{sample}/bsj.bed",
        genome = config["ref"]["genome"],
        ann    = "01_data/annotation/{}.txt".format(config["ref"]["prefix"])
    output: "05_output/circexplorer/{aligner}/{sample}/circRNA_known.txt"
    log: "logs/circexplorer/{aligner}/{sample}_annotate.log"
    conda: "../envs/circ.yaml"
    shell :
        """
           CIRCexplorer2 annotate -r {input.ann} -g {input.genome} -b {input.bsj} -o {output} > {log}
        """

rule gtf2genepred:
    input : 
        gtf = config["ref"]["annotation"]
    output :
        pred = "01_data/annotation/{}.genepred".format(config["ref"]["prefix"]),
        txt = "01_data/annotation/{}.txt".format(config["ref"]["prefix"])
    conda : "../envs/circ.yaml"
    shell : 
        """
            gtfToGenePred -geneNameAsName2 -genePredExt {input.gtf} {output.pred};
            awk '{{printf $12"\t"}}{{for(i=1;i<11;i++){{if(i==10){{printf $i}} else{{printf $i"\\t"}}}}; {{printf "\\n"}}}}'  {output.pred} > {output.txt}
        """
