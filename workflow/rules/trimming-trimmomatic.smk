# =================================================================================================
#     Trimming
# =================================================================================================


rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        (
            "trimming/{sample}-{unit}.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.fastq.gz")
        ),
        # trimlog="trimming/{sample}-{unit}-se.trimlog.log"
        touch("trimming/{sample}-{unit}-se.done"),
    params:
        # extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        extra=config["params"]["trimmomatic"]["se"]["extra"],
        trimmer=config["params"]["trimmomatic"]["se"]["trimmer"],
    threads: config["params"]["trimmomatic"]["threads"]
    log:
        "logs/trimming/trimmomatic/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/trimmomatic/{sample}-{unit}.log"
    conda:
        # yet another missing dependency in the original wrapper...
        "../envs/trimmomatic.yaml"
    wrapper:
        "v3.13.6/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=(
            "trimming/{sample}-{unit}.1.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.1.fastq.gz")
        ),
        r2=(
            "trimming/{sample}-{unit}.2.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.2.fastq.gz")
        ),
        r1_unpaired=(
            "trimming/{sample}-{unit}.1.unpaired.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.1.unpaired.fastq.gz")
        ),
        r2_unpaired=(
            "trimming/{sample}-{unit}.2.unpaired.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.2.unpaired.fastq.gz")
        ),
        # trimlog="trimming/{sample}-{unit}-pe.trimlog.log"
        done=touch("trimming/{sample}-{unit}-pe.done"),
    params:
        # extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        extra=config["params"]["trimmomatic"]["se"]["extra"],
        trimmer=config["params"]["trimmomatic"]["pe"]["trimmer"],
    threads: config["params"]["trimmomatic"]["threads"]
    log:
        "logs/trimming/trimmomatic/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/trimmomatic/{sample}-{unit}.log"
    conda:
        # yet another missing dependency in the original wrapper...
        "../envs/trimmomatic.yaml"
    wrapper:
        "v3.13.6/bio/trimmomatic/pe"


# =================================================================================================
#     Trimming Results
# =================================================================================================


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if is_single_end(wildcards.sample, wildcards.unit):
        # single end sample
        return [
            "trimming/{sample}-{unit}.fastq.gz".format(sample=wildcards.sample, unit=wildcards.unit)
        ]
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        raise Exception(
            "Trimming tool 'trimmomatic' cannot be used with the option 'merge-paired-end-reads'"
        )
    else:
        # paired-end sample
        return expand(
            "trimming/{sample}-{unit}.{pair}.fastq.gz",
            pair=[1, 2],
            sample=wildcards.sample,
            unit=wildcards.unit,
        )


# MultiQC expects the normal stdout log files from trimmomatic, but as we use a wrapper for the latter,
# we cannot also declare the log files as output files, because snakemake...
# Instead, we copy them afterwards. This is dirty, but that's how the snake rolls...
rule trimmomatic_multiqc_log:
    input:
        # Take the trimming result as dummy input, so that this rule here is always executed afterwards
        get_trimmed_reads,
    output:
        "trimming/{sample}-{unit}.trimlog.log",
    shell:
        "cp logs/trimming/trimmomatic/{wildcards.sample}-{wildcards.unit}.log {output}"


localrules:
    trimmomatic_multiqc_log,


def get_trimming_report(sample, unit):
    """Get the report needed for MultiQC."""
    if is_single_end(sample, unit):
        # single end sample
        return "trimming/" + sample + "-" + unit + ".trimlog.log"
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        raise Exception(
            "Trimming tool 'trimmomatic' cannot be used with the option 'merge-paired-end-reads'"
        )
    else:
        # paired-end sample
        return "trimming/" + sample + "-" + unit + ".trimlog.log"
