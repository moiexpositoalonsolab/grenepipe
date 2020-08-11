# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        r1=temp("trimmed/{sample}-{unit}-trimmed.fastq.gz"),
        settings="trimmed/{sample}-{unit}-trimmed-se.settings"
    params:
        extra="--gzip",
        params=config["params"]["adapterremoval"]["se"],
        basename="trimmed/{sample}-{unit}-trimmed"
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
        r1=temp("trimmed/{sample}-{unit}-trimmed-pair1.fastq.gz"),
        r2=temp("trimmed/{sample}-{unit}-trimmed-pair2.fastq.gz"),
        settings="trimmed/{sample}-{unit}-trimmed-pe.settings"
    params:
        extra="--gzip",
        params=config["params"]["adapterremoval"]["pe"],
        basename="trimmed/{sample}-{unit}-trimmed"
    threads:
        config["params"]["adapterremoval"]["threads"]
    log:
        "logs/adapterremoval/{sample}-{unit}.log"
    benchmark:
        "benchmarks/adapterremoval/{sample}-{unit}.bench.log"
    conda:
        "../envs/adapterremoval.yaml"
    shell:
        "AdapterRemoval --file1 {input.r1} --file2 {input.r2} {params.extra} {params.params} --threads {threads} "
        "--basename {params.basename} --output1 {output.r1} --output2 {output.r2} --settings {output.settings} > {log} 2>&1"

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
        return [ "trimmed/{sample}-{unit}-trimmed.fastq.gz".format(**wildcards) ]
    else:
        # paired-end sample
        return expand("trimmed/{sample}-{unit}-trimmed-pair{group}.fastq.gz", group=[1, 2], **wildcards)
