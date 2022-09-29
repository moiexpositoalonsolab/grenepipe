# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        merged=(
            "trimmed/{sample}-{unit}-merged.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}-merged.fastq.gz")
        ),
        fastq1=(
            "trimmed/{sample}-{unit}.1.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.1.fastq.gz")
        ),
        fastq2=(
            "trimmed/{sample}-{unit}.2.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.2.fastq.gz")
        ),
        fastq1_discarded=(
            "trimmed/{sample}-{unit}.1.discarded.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.1.discarded.fastq.gz")
        ),
        fastq2_discarded=(
            "trimmed/{sample}-{unit}.2.discarded.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.2.discarded.fastq.gz")
        ),
        done=touch("trimmed/{sample}-{unit}.done")
    params:
        extra = config["params"]["seqprep"]["extra"]
    log:
        "logs/seqprep/{sample}-{unit}.log"
    benchmark:
        "benchmarks/seqprep/{sample}-{unit}.bench.log"
    conda:
        "../envs/seqprep.yaml"
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
        raise Exception(
            "Trimming tool 'seqprep' cannot be used for single-end reads"
        )
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return [ "trimmed/{sample}-{unit}-merged.fastq.gz".format(
            sample=wildcards.sample, unit=wildcards.unit
        )]
    else:
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{pair}.fastq.gz",
            pair=[1, 2], sample=wildcards.sample, unit=wildcards.unit
        )

# MultiQC does not support SeqPrep. Empty output.
def get_trimming_report(sample, unit):
    """Get the report needed for MultiQC."""
    if is_single_end(sample, unit):
        # single end sample
        raise Exception(
            "Trimming tool 'seqprep' cannot be used for single-end reads"
        )
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return []
    else:
        # paired-end sample
        return []
