# =================================================================================================
#     Genome Indexing
# =================================================================================================

# Bowtie needs its own index of the reference genome. Of course it does.
rule bowtie2_index:
    input:
        config["data"]["reference"]["genome"]
    output:
        output1=config["data"]["reference"]["genome"] + ".1.bt2",
        output2=config["data"]["reference"]["genome"] + ".2.bt2",
        output3=config["data"]["reference"]["genome"] + ".3.bt2",
        output4=config["data"]["reference"]["genome"] + ".4.bt2",
        outputrev1=config["data"]["reference"]["genome"] + ".rev.1.bt2",
        outputrev2=config["data"]["reference"]["genome"] + ".rev.2.bt2"
    params:
        # Bowtie expects the prefix, and creates the above output files automatically.
        # So, let's do some snakemake magic to make this work.
        # Let's hope that they do not change the naming of their ouput files in the future.
        basename=config["data"]["reference"]["genome"]
    conda:
        "../envs/bowtie2.yaml"
    log:
        os.path.dirname(config["data"]["reference"]["genome"]) + "/logs/bowtie2_index.log"
    shell:
        "bowtie2-build {input} {params.basename} > {log} 2>&1"

# =================================================================================================
#     Read Mapping
# =================================================================================================

rule map_reads:
    input:
        sample=get_trimmed_reads,
        bowtie_idx1=config["data"]["reference"]["genome"] + ".1.bt2",
        bowtie_idx2=config["data"]["reference"]["genome"] + ".2.bt2",
        bowtie_idx3=config["data"]["reference"]["genome"] + ".3.bt2",
        bowtie_idx4=config["data"]["reference"]["genome"] + ".4.bt2",
        bowtie_rev1=config["data"]["reference"]["genome"] + ".rev.1.bt2",
        bowtie_rev2=config["data"]["reference"]["genome"] + ".rev.2.bt2"
    output:
        pipe(config["rundir"] + "mapped/{sample}-{unit}.bam")
    params:
        # Prefix of reference genome index (built with bowtie2-build)
        index=config["data"]["reference"]["genome"],

        # Optional parameters.
        # We need the read group tags, including `ID` and `SM`, as downstream tools use these.
        extra="--rg-id {sample} --rg SM:{sample}"
    threads:
        # Use at least two threads
        config["params"]["bowtie2"]["threads"]
    log:
        config["rundir"] + "logs/bowtie2/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/bowtie2/{sample}-{unit}.bench.log"
    wrapper:
        "0.58.0/bio/bowtie2/align"

# The bowtie wrapper above is different from the bwa mem wrapper, and does not sort.
# Why is everything so inconsistent? Anyway, that means we have to do the sorting step here.
# At least, we can pipe the files from above to here, so this should not slow us down.
rule sort_reads:
    input:
        config["rundir"] + "mapped/{sample}-{unit}.bam"
    output:
        config["rundir"] + "mapped/{sample}-{unit}.sorted.bam"
    params:
        "-m 4G"
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@. Weird flex, but okay.
    log:
        config["rundir"] + "logs/samtools/sort/{sample}-{unit}.log"
    wrapper:
        "0.58.0/bio/samtools/sort"
