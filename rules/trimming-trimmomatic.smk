# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        "trimmed/{sample}-{unit}.fastq.gz",
        trimlog="trimmed/{sample}-{unit}-se.trimlog.log"
    params:
        # extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["se"]
    threads:
        config["params"]["trimmomatic"]["threads"]
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    benchmark:
        "benchmarks/trimmomatic/{sample}-{unit}.bench.log"
    wrapper:
        "0.51.3/bio/trimmomatic/se"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1="trimmed/{sample}-{unit}.1.fastq.gz",
        r2="trimmed/{sample}-{unit}.2.fastq.gz",
        r1_unpaired="trimmed/{sample}-{unit}.1.unpaired.fastq.gz",
        r2_unpaired="trimmed/{sample}-{unit}.2.unpaired.fastq.gz",
        trimlog="trimmed/{sample}-{unit}-pe.trimlog.log"
    params:
        # extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    threads:
        config["params"]["trimmomatic"]["threads"]
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    benchmark:
        "benchmarks/trimmomatic/{sample}-{unit}.bench.log"
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
        return [ "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards) ]
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        raise Exception("Trimming tool 'trimmomatic' cannot be used with the option 'merge-paired-end-reads'")
    else:
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards)

# MultiQC expects the normal stdout log files from trimmomatic, but as we use a wrapper for the latter,
# we cannot also declare the log files as output files, because snakemake...
# Instead, we copy them afterwards. This is dirty, but that's how the snake rolls...
rule trimmomatic_multiqc_log:
    input:
        get_trimmed_reads
    output:
        "trimmed/{sample}-{unit}.trimlog.log"
    shell:
        "cp logs/trimmomatic/{wildcards.sample}-{wildcards.unit}.log {output}"

def get_trimming_report(sample, unit):
    """Get the report needed for MultiQC."""
    if is_single_end(sample, unit):
        # single end sample
        return "trimmed/" + sample + "-" + unit + "-se.trimlog.log"
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        raise Exception("Trimming tool 'trimmomatic' cannot be used with the option 'merge-paired-end-reads'")
    else:
        # paired-end sample
        return "trimmed/" + sample + "-" + unit + "-pe.trimlog.log"
