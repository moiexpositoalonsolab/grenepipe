# =================================================================================================
#     Read Mapping
# =================================================================================================

# Apprently, samtools does not create the tmp dir correclty, so we need to take care of this...
if len(config["params"]["samtools"]["temp-dir"]) > 0:
    os.makedirs(config["params"]["samtools"]["temp-dir"], exist_ok=True)


def get_bwa_mem_extra(wildcards):
    rg_tags = "\\t".join(get_read_group_tags(wildcards))
    extra = "-R '@RG\\t" + rg_tags + "' " + config["params"]["bwamem"]["extra"]
    return extra


rule map_reads:
    input:
        reads=get_trimmed_reads,
        ref=config["data"]["reference-genome"],
        idx=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
    output:
        (
            "mapping/sorted/{sample}-{unit}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapping/sorted/{sample}-{unit}.bam")
        ),
        touch("mapping/sorted/{sample}-{unit}.done"),
    params:
        index=config["data"]["reference-genome"],
        extra=get_bwa_mem_extra,
        # Sort as we need it. Samtools provided via the two params for new and old wrapper versions.
        sort="samtools",
        sorting="samtools",
        sort_order="coordinate",
        sort_extra=config["params"]["samtools"]["sort"],
        tmp_dir=config["params"]["samtools"]["temp-dir"],
    group:
        "mapping"
    log:
        "logs/mapping/bwa-mem/{sample}-{unit}.log",
    benchmark:
        "benchmarks/mapping/bwa-mem/{sample}-{unit}.log"
    threads: config["params"]["bwamem"]["threads"]
    conda:
        "../envs/bwa.yaml"
    # resources:
    # Increase time limit in factors of 2h, if the job fails due to time limit.
    # time = lambda wildcards, input, threads, attempt: int(120 * int(attempt))

    # This wrapper version uses a proper tmp dir, so that the below shadow rule is not needed.
    # It caused trouble when running large cluster jobs with high number of parallel jobs,
    # as the number of symlinks created for the shadow directory crashed our cluster max file
    # count limit...
    wrapper:
        "0.80.0/bio/bwa/mem"


# samtools sort creates tmp files that are not cleaned up when the cluster job runs out
# of time, but which cause samtools to immediately terminate if called again, meaning
# that we cannot run it again with more time. We hence use a full shadow directory.
# We experimented with all kinds of other solutions, such as setting the `-T` option of
# `samtools sort` via the `sort_extra` param, and creating and deleting tmp directories for that,
# but that fails due to issues with snakemake handling directories instead of files...
# The snakemake shadow rule here gets the job done.
# shadow:
#     "full"
# wrapper:
#     "0.51.3/bio/bwa/mem"
