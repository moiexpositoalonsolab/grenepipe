# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        r1="trimmed/{sample}-{unit}.fastq.gz",
        settings="trimmed/{sample}-{unit}.se.settings"
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
        r1="trimmed/{sample}-{unit}.pair1.fastq.gz",
        r2="trimmed/{sample}-{unit}.pair2.fastq.gz",
        settings="trimmed/{sample}-{unit}.pe.settings"
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
        pair1_truncated     = "trimmed/{sample}-{unit}.pair1.truncated.gz",
        pair2_truncated     = "trimmed/{sample}-{unit}.pair2.truncated.gz",
        singleton_truncated = "trimmed/{sample}-{unit}.singleton.truncated.gz",
        collapsed           = "trimmed/{sample}-{unit}.collapsed.gz",
        collapsed_truncated = "trimmed/{sample}-{unit}.collapsed.truncated.gz",
        discarded           = "trimmed/{sample}-{unit}.discarded.gz",
        settings            = "trimmed/{sample}-{unit}.pe-merged.settings"
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
        return [ "trimmed/{sample}-{unit}.collapsed.gz".format(**wildcards) ]
    else:
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.pair{group}.fastq.gz", group=[1, 2], **wildcards)

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
