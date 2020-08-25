# =================================================================================================
#     Trimming
# =================================================================================================

# The fastp wrapper is different from the other trimming wrappers in the the snakemake wrapper
# respository, because consistency is just not their strength... So, we have to provide an extra
# function to unpack the fastq file names into a list for us.
def unpack_fastp_files(wildcards):
    return list(get_fastq(wildcards).values())

rule trim_reads_se:
    input:
        sample=unpack_fastp_files
    output:
        trimmed=temp("trimmed/{sample}-{unit}.fastq.gz"),
        html="trimmed/{sample}-{unit}-se-fastp.html",
        json="trimmed/{sample}-{unit}-se-fastp.json"
    log:
        "logs/fastp/{sample}-{unit}.log"
    benchmark:
        "benchmarks/fastp/{sample}-{unit}.bench.log"
    params:
        extra=config["params"]["fastp"]["se"]
    threads:
        config["params"]["fastp"]["threads"]
    wrapper:
        "0.64.0/bio/fastp"

rule trim_reads_pe:
    input:
        sample=unpack_fastp_files
    output:
        trimmed=temp(["trimmed/{sample}-{unit}.1.fastq.gz", "trimmed/{sample}-{unit}.2.fastq.gz"]),
        html="trimmed/{sample}-{unit}-pe-fastp.html",
        json="trimmed/{sample}-{unit}-pe-fastp.json"
    log:
        "logs/fastp/{sample}-{unit}.log"
    benchmark:
        "benchmarks/fastp/{sample}-{unit}.bench.log"
    params:
        extra=config["params"]["fastp"]["pe"]
    threads:
        config["params"]["fastp"]["threads"]
    wrapper:
        "0.64.0/bio/fastp"

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
