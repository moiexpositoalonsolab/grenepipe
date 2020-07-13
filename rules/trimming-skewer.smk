# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp("trimmed/{sample}-{unit}-trimmed.fastq.gz")
    params:
        extra="--format sanger --compress",
        params=config["params"]["skewer"]["se"],
        outpref="trimmed/{sample}-{unit}"
    threads:
        config["params"]["skewer"]["threads"]
    log:
        "logs/skewer/{sample}-{unit}.log"
    benchmark:
        "benchmarks/skewer/{sample}-{unit}.bench.log"
    conda:
        "../envs/skewer.yaml"
    shell:
        "skewer {params.extra} {params.params} --threads {threads} --output {params.outpref} "
        "{input.r1} > {log} 2>&1"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp("trimmed/{sample}-{unit}-trimmed-pair1.fastq.gz"),
        r2=temp("trimmed/{sample}-{unit}-trimmed-pair2.fastq.gz")
    params:
        extra="--format sanger --compress",
        params=config["params"]["skewer"]["pe"],
        outpref="trimmed/{sample}-{unit}"
    threads:
        config["params"]["skewer"]["threads"]
    log:
        "logs/skewer/{sample}-{unit}.log"
    benchmark:
        "benchmarks/skewer/{sample}-{unit}.bench.log"
    conda:
        "../envs/skewer.yaml"
    shell:
        "skewer {params.extra} {params.params} --threads {threads} --output {params.outpref} "
        "{input.r1} {input.r2} > {log} 2>&1"

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
