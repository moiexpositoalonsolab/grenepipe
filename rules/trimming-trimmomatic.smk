# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp(config["rundir"] + "trimmed/{sample}-{unit}.fastq.gz")
    params:
        extra="",
        **config["params"]["trimmomatic"]["se"]
    threads:
        config["params"]["trimmomatic"]["threads"]
    log:
        config["rundir"] + "logs/trimmomatic/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/trimmomatic/{sample}-{unit}.bench.log"
    wrapper:
        "0.51.3/bio/trimmomatic/se"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp(config["rundir"] + "trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp(config["rundir"] + "trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp(config["rundir"] + "trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp(config["rundir"] + "trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog=config["rundir"] + "trimmed/{sample}-{unit}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    threads:
        config["params"]["trimmomatic"]["threads"]
    log:
        config["rundir"] + "logs/trimmomatic/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/trimmomatic/{sample}-{unit}.bench.log"
    wrapper:
        "0.51.3/bio/trimmomatic/pe"

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
        return [ config["rundir"] + "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards) ]
    else:
        # paired-end sample
        return expand(config["rundir"] + "trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards)
