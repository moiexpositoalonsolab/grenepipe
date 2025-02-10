# =================================================================================================
#     Trimming
# =================================================================================================


rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        fastq=(
            "trimming/{sample}-{unit}.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.fastq.gz")
        ),
        qc="trimming/{sample}-{unit}.qc-se.txt",
        done=touch("trimming/{sample}-{unit}.fastq.gz.done"),
    params:
        adapters=config["params"]["cutadapt"]["se"]["adapters"],
        extra=config["params"]["cutadapt"]["se"]["extra"],
    threads: config["params"]["cutadapt"]["threads"]
    log:
        "logs/trimming/cutadapt/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/cutadapt/{sample}-{unit}.log"
    conda:
        # yet another missing dependency in the original wrapper...
        "../envs/cutadapt.yaml"
    wrapper:
        "0.74.0/bio/cutadapt/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        fastq1=(
            "trimming/{sample}-{unit}.1.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.1.fastq.gz")
        ),
        fastq2=(
            "trimming/{sample}-{unit}.2.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.2.fastq.gz")
        ),
        qc="trimming/{sample}-{unit}.qc-pe.txt",
        done1=touch("trimming/{sample}-{unit}.1.fastq.gz.done"),
        done2=touch("trimming/{sample}-{unit}.2.fastq.gz.done"),
    params:
        adapters=config["params"]["cutadapt"]["pe"]["adapters"],
        extra=config["params"]["cutadapt"]["pe"]["extra"],
    threads: config["params"]["cutadapt"]["threads"]
    log:
        "logs/trimming/cutadapt/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/cutadapt/{sample}-{unit}.log"
    conda:
        # yet another missing dependency in the original wrapper...
        "../envs/cutadapt.yaml"
    wrapper:
        "0.74.0/bio/cutadapt/pe"


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
            "Trimming tool 'cutadapt' cannot be used with the option 'merge-paired-end-reads'"
        )
    else:
        # paired-end sample
        return expand(
            "trimming/{sample}-{unit}.{pair}.fastq.gz",
            sample=wildcards.sample,
            unit=wildcards.unit,
            pair=[1, 2],
        )


def get_trimmed_reads_done(wildcards):
    files = get_trimmed_reads(wildcards)
    return [f + ".done" for f in files]


def get_trimming_report(sample, unit):
    """Get the report needed for MultiQC."""
    if is_single_end(sample, unit):
        # single end sample
        return "trimming/" + sample + "-" + unit + ".qc-se.txt"
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        raise Exception(
            "Trimming tool 'cutadapt' cannot be used with the option 'merge-paired-end-reads'"
        )
    else:
        # paired-end sample
        return "trimming/" + sample + "-" + unit + ".qc-pe.txt"
