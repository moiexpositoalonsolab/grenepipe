# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        fastq=temp("trimmed/{sample}-{unit}.fastq.gz"),
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        config["params"]["cutadapt"]["se"]
    threads:
        config["params"]["cutadapt"]["threads"]
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    benchmark:
        "benchmarks/cutadapt/{sample}-{unit}.bench.log"
    wrapper:
        "0.63.0/bio/cutadapt/se"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        fastq1=temp("trimmed/{sample}-{unit}.1.fastq.gz"),
        fastq2=temp("trimmed/{sample}-{unit}.2.fastq.gz"),
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        adapters = config["params"]["cutadapt"]["pe"]["adapters"],
        others = config["params"]["cutadapt"]["pe"]["others"]
    threads:
        config["params"]["cutadapt"]["threads"]
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    benchmark:
        "benchmarks/cutadapt/{sample}-{unit}.bench.log"
    wrapper:
        "0.63.0/bio/cutadapt/pe"

# =================================================================================================
#     Trimming Results
# =================================================================================================

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if is_single_end(**wildcards):
        # single end sample
        return [ "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards) ]
    else:
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards)
