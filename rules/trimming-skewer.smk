# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp(config["rundir"] + "trimmed/{sample}-{unit}-trimmed.fastq.gz")
    params:
        extra="--format sanger --compress",
        params=config["params"]["skewer"]["se"],
        outpref=config["rundir"] + "trimmed/{sample}-{unit}"
    threads:
        config["params"]["skewer"]["threads"]
    log:
        config["rundir"] + "logs/skewer/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/skewer/{sample}-{unit}.bench.log"
    conda:
        "../envs/skewer.yaml"
    shell:
        "skewer {params.extra} {params.params} --threads {threads} --output {params.outpref} "
        "{input.r1} > {log} 2>&1"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp(config["rundir"] + "trimmed/{sample}-{unit}-trimmed-pair1.fastq.gz"),
        r2=temp(config["rundir"] + "trimmed/{sample}-{unit}-trimmed-pair2.fastq.gz")
    params:
        extra="--format sanger --compress",
        params=config["params"]["skewer"]["pe"],
        outpref=config["rundir"] + "trimmed/{sample}-{unit}"
    threads:
        config["params"]["skewer"]["threads"]
    log:
        config["rundir"] + "logs/skewer/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/skewer/{sample}-{unit}.bench.log"
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
        return [ config["rundir"] + "trimmed/{sample}-{unit}-trimmed.fastq.gz".format(**wildcards) ]
    else:
        # paired-end sample
        return expand(config["rundir"] + "trimmed/{sample}-{unit}-trimmed-pair{group}.fastq.gz", group=[1, 2], **wildcards)
