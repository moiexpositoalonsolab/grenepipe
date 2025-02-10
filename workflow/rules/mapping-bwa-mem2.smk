# =================================================================================================
#     Read Mapping
# =================================================================================================

# Apprently, samtools does not create the tmp dir correclty, so we need to take care of this...
if len(config["params"]["samtools"]["temp-dir"]) > 0:
    os.makedirs(config["params"]["samtools"]["temp-dir"], exist_ok=True)


def get_bwa_mem2_extra(wildcards):
    rg_tags = "\\t".join(get_read_group_tags(wildcards))
    extra = "-R '@RG\\t" + rg_tags + "' " + config["params"]["bwamem2"]["extra"]
    return extra


rule map_reads:
    input:
        reads=get_trimmed_reads,
        reads_done=get_trimmed_reads_done,
        ref=config["data"]["reference-genome"],
        # Somehow, the wrapper expects the index extensions to be given,
        # instead of the underlying fasta file... Well, so let's do that.
        # We provide the fasta above as well; it's not used,
        # but might be important as a rule dependency so that it is present.
        idx=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["0123", "amb", "ann", "bwt.2bit.64", "pac"],
        ),
    output:
        (
            "mapping/sorted/{sample}-{unit}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapping/sorted/{sample}-{unit}.bam")
        ),
        touch("mapping/sorted/{sample}-{unit}.bam.done"),
    params:
        extra=get_bwa_mem2_extra,
        # Sort as we need it.
        sort="samtools",
        sort_order="coordinate",
        sort_extra=config["params"]["samtools"]["sort"],
        tmp_dir=config["params"]["samtools"]["temp-dir"],
    group:
        "mapping"
    log:
        "logs/mapping/bwa-mem2/{sample}-{unit}.log",
    benchmark:
        "benchmarks/mapping/bwa-mem2/{sample}-{unit}.log"
    threads: config["params"]["bwamem2"]["threads"]
    conda:
        # As always, we need our own env here that overwrites the python/pandas/numpy stack
        # to make sure that we do not run into a version conflict.
        "../envs/bwa-mem2.yaml"
    script:
        # We use our own version of the wrapper here, as that wrapper misses temp dirs for
        # samtools sort, causing all kinds of trouble...
        "../scripts/bwa-mem2-mem.py"


# wrapper:
#     "0.78.0/bio/bwa-mem2/mem"
