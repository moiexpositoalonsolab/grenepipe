# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        r1=(
            "trimmed/{sample}-{unit}.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.fastq.gz")
        ),
        settings="trimmed/{sample}-{unit}.se.settings",
        done=touch("trimmed/{sample}-{unit}.se.done")
    params:
        extra="--gzip",
        params=config["params"]["adapterremoval"]["se"],
        basename="trimmed/{sample}-{unit}"
    threads:
        config["params"]["adapterremoval"]["threads"]
    log:
        "logs/adapterremoval/{sample}-{unit}.log"
    benchmark:
        "benchmarks/adapterremoval/{sample}-{unit}.bench.log"
    conda:
        "../envs/adapterremoval.yaml"
    shell:
        "AdapterRemoval --file1 {input.r1} {params.extra} {params.params} --threads {threads} "
        "--basename {params.basename} --output1 {output.r1} --settings {output.settings} > {log} 2>&1"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=(
            "trimmed/{sample}-{unit}.pair1.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.pair1.fastq.gz")
        ),
        r2=(
            "trimmed/{sample}-{unit}.pair2.fastq.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.pair2.fastq.gz")
        ),
        settings="trimmed/{sample}-{unit}.pe.settings",
        done=touch("trimmed/{sample}-{unit}.pe.done")
    params:
        extra="--gzip",
        params=config["params"]["adapterremoval"]["pe"],
        basename="trimmed/{sample}-{unit}"
    threads:
        config["params"]["adapterremoval"]["threads"]
    log:
        "logs/adapterremoval/{sample}-{unit}.log"
    benchmark:
        "benchmarks/adapterremoval/{sample}-{unit}.bench.log"
    conda:
        "../envs/adapterremoval.yaml"
    shell:
        "AdapterRemoval --file1 {input.r1} --file2 {input.r2} {params.extra} {params.params} "
        "--threads {threads} --basename {params.basename} --output1 {output.r1} --output2 {output.r2} "
        "--settings {output.settings} > {log} 2>&1"

rule trim_reads_pe_merged:
    input:
        unpack(get_fastq)
    output:
        # We mostly use the default output file names of AdapterRemoval here,
        # except for the settings file, which we have to rename so that its name is distinct from
        # the settings files of the other two rules above.
        pair1_truncated = (
            "trimmed/{sample}-{unit}.pair1.truncated.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.pair1.truncated.gz")
        ),
        pair2_truncated = (
            "trimmed/{sample}-{unit}.pair2.truncated.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.pair2.truncated.gz")
        ),
        singleton_truncated = (
            "trimmed/{sample}-{unit}.singleton.truncated.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.singleton.truncated.gz")
        ),
        collapsed = (
            "trimmed/{sample}-{unit}.collapsed.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.collapsed.gz")
        ),
        collapsed_truncated = (
            "trimmed/{sample}-{unit}.collapsed.truncated.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.collapsed.truncated.gz")
        ),
        discarded = (
            "trimmed/{sample}-{unit}.discarded.gz"
            if config["settings"]["keep-intermediate"]["trimming"]
            else temp("trimmed/{sample}-{unit}.discarded.gz")
        ),
        settings = "trimmed/{sample}-{unit}.pe-merged.settings",
        done = touch("trimmed/{sample}-{unit}.pe-merged.done")
    params:
        extra="--gzip --collapse",
        params=config["params"]["adapterremoval"]["pe"],
        basename="trimmed/{sample}-{unit}"
    threads:
        config["params"]["adapterremoval"]["threads"]
    log:
        "logs/adapterremoval/{sample}-{unit}.log"
    benchmark:
        "benchmarks/adapterremoval/{sample}-{unit}.bench.log"
    conda:
        "../envs/adapterremoval.yaml"
    shell:
        "AdapterRemoval --file1 {input.r1} --file2 {input.r2} {params.extra} {params.params} "
        "--threads {threads} --basename {params.basename} --settings {output.settings} > {log} 2>&1"

# =================================================================================================
#     Trimming Results
# =================================================================================================

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if is_single_end(wildcards.sample, wildcards.unit):
        # single end sample
        return [ "trimmed/{sample}-{unit}.fastq.gz".format(
            sample=wildcards.sample, unit=wildcards.unit
        )]
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return [ "trimmed/{sample}-{unit}.collapsed.gz".format(
            sample=wildcards.sample, unit=wildcards.unit
        )]
    else:
        # paired-end sample
        return expand(
            "trimmed/{sample}-{unit}.pair{pair}.fastq.gz",
            pair=[1, 2], sample=wildcards.sample, unit=wildcards.unit
        )

def get_trimming_report(sample, unit):
    """Get the report needed for MultiQC."""
    if is_single_end(sample, unit):
        # single end sample
        return "trimmed/" + sample + "-" + unit + ".se.settings"
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return "trimmed/" + sample + "-" + unit + ".pe-merged.settings"
    else:
        # paired-end sample
        return "trimmed/" + sample + "-" + unit + ".pe.settings"
