import platform

# =================================================================================================
#     Trimming
# =================================================================================================

if platform.system() == "Darwin":
    raise Exception(
        "Trimming tool seqprep is not available for MacOS. Please use a different trimming tool."
    )


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        merged=(
            "trimming/{sample}-{unit}-merged.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}-merged.fastq.gz")
        ),
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
        fastq1_discarded=(
            "trimming/{sample}-{unit}.1.discarded.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.1.discarded.fastq.gz")
        ),
        fastq2_discarded=(
            "trimming/{sample}-{unit}.2.discarded.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimming/{sample}-{unit}.2.discarded.fastq.gz")
        ),
        done=touch("trimming/{sample}-{unit}.done"),
    params:
        extra=config["params"]["seqprep"]["extra"],
    log:
        "logs/trimming/seqprep/{sample}-{unit}.log",
    benchmark:
        "benchmarks/trimming/seqprep/{sample}-{unit}.log"
    conda:
        "../envs/seqprep-linux.yaml"
    shell:
        "SeqPrep "
        "-f {input.r1} -r {input.r2} "
        "-1 {output.fastq1} -2 {output.fastq2} "
        "-3 {output.fastq1_discarded} -4 {output.fastq2_discarded} "
        "-s {output.merged} "
        "{params.extra} > {log} 2>&1"


# =================================================================================================
#     Trimming Results
# =================================================================================================


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if is_single_end(wildcards.sample, wildcards.unit):
        # single end sample
        raise Exception("Trimming tool 'seqprep' cannot be used for single-end reads")
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return [
            "trimming/{sample}-{unit}-merged.fastq.gz".format(
                sample=wildcards.sample, unit=wildcards.unit
            )
        ]
    else:
        # paired-end sample
        return expand(
            "trimming/{sample}-{unit}.{pair}.fastq.gz",
            pair=[1, 2],
            sample=wildcards.sample,
            unit=wildcards.unit,
        )


# MultiQC does not support SeqPrep. Empty output.
def get_trimming_report(sample, unit):
    """Get the report needed for MultiQC."""
    if is_single_end(sample, unit):
        # single end sample
        raise Exception("Trimming tool 'seqprep' cannot be used for single-end reads")
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return []
    else:
        # paired-end sample
        return []
