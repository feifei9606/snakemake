def is_single_end(sample):
    return pd.isnull(units[units['sample'] == sample]['fq2'])[0]

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if peFlag:
        # paired-end sample
        return expand("03_trim/{sample}.{group}.trim.fastq.gz",
                      group=["R1", "R2"], **wildcards)
    # single end sample
    return expand("03_trim/{sample}.trim.fastq.gz".format(**wildcards))