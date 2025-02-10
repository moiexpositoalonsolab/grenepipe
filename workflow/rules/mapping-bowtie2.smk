# =================================================================================================
#     Genome Indexing
# =================================================================================================


# Bowtie needs its own index of the reference genome. Of course it does.
rule bowtie2_index:
    input:
        config["data"]["reference-genome"],
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
    output:
        expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"],
        ),
        done=touch(config["data"]["reference-genome"] + ".done")
    params:
        # Bowtie expects the prefix, and creates the above output files automatically.
        # So, let's do some snakemake magic to make this work.
        # Let's hope that they do not change the naming of their ouput files in the future.
        basename=config["data"]["reference-genome"],
    conda:
        "../envs/bowtie2.yaml"
    log:
        os.path.dirname(config["data"]["reference-genome"]) + "/logs/bowtie2_index.log",
    shell:
        "bowtie2-build {input[0]} {params.basename} > {log} 2>&1"


# =================================================================================================
#     Read Mapping
# =================================================================================================

# Apprently, samtools does not create the tmp dir correclty, so we need to take care of this...
if len(config["params"]["samtools"]["temp-dir"]) > 0:
    os.makedirs(config["params"]["samtools"]["temp-dir"], exist_ok=True)


def get_bowtie2_extra(wildcards):
    extra = r"--rg-id \"" + wildcards.sample + r"\""
    rg_tags = get_read_group_tags(wildcards)
    for tag in rg_tags:
        extra += r" --rg \"" + tag + r"\""
    extra += " " + config["params"]["bowtie2"]["extra"]
    return extra


rule map_reads:
    input:
        sample=get_trimmed_reads,
        sample_done=get_trimmed_reads_done,
        indices=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2", "done"],
        ),
    output:
        pipe("mapping/mapped/{sample}-{unit}.bam"),
        # touch("mapping/mapped/{sample}-{unit}.bam.done"),
    params:
        # Prefix of reference genome index (built with bowtie2-build above)
        index=config["data"]["reference-genome"],
        extra=get_bowtie2_extra,
    # Use at least two threads
    threads: config["params"]["bowtie2"]["threads"]
    # resources:
    # Increase time limit in factors of 2h, if the job fails due to time limit.
    # time = lambda wildcards, input, threads, attempt: int(120 * int(attempt))
    log:
        "logs/mapping/bowtie2/{sample}-{unit}.log",
    benchmark:
        "benchmarks/mapping/bowtie2/{sample}-{unit}.log"
    group:
        "mapping"
    conda:
        # We need our own conda environment here, as the wrapper environment is not specifying
        # the exact version of libtbb, which in its latest version does not work with bowtie2...
        "../envs/bowtie2.yaml"
    wrapper:
        "0.74.0/bio/bowtie2/align"


# The bowtie wrapper above is different from the bwa mem wrapper, and does not sort.
# Why is everything so inconsistent? Anyway, that means we have to do the sorting step here.
# At least, we can pipe the files from above to here, so this should not slow us down.
rule sort_reads:
    input:
        "mapping/mapped/{sample}-{unit}.bam",
        "mapping/mapped/{sample}-{unit}.bam.done",
    output:
        (
            "mapping/sorted/{sample}-{unit}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapping/sorted/{sample}-{unit}.bam")
        ),
        touch("mapping/sorted/{sample}-{unit}.bam.done"),
    params:
        extra=config["params"]["samtools"]["sort"],
        tmp_dir=config["params"]["samtools"]["temp-dir"],
    # Samtools takes additional threads through its option -@
    threads: 1  # This value - 1 will be sent to -@. Weird flex, but okay.
    log:
        "logs/mapping/samtools-sort/{sample}-{unit}.log",
    group:
        "mapping"
    conda:
        "../envs/samtools.yaml"
    wrapper:
        "0.80.0/bio/samtools/sort"
