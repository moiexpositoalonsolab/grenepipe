# =================================================================================================
#     Merge bams
# =================================================================================================

rule mpileup_merge_all:
    input:
        get_all_bams()
    output:
        pipe("mpileup/all.merged.bam"),
        touch("mpileup/all.merged.done")
    params:
        # Need file overwrite flag, see above.
        extra=config["params"]["samtools"]["merge"] + " -f"
    threads:
        # Samtools takes additional threads through its option -@
        # This value - 1 will be sent to -@
        config["params"]["samtools"]["merge-threads"]
    log:
        "logs/samtools/mpileup/merge-all.log"
    wrapper:
        "v3.13.6/bio/samtools/merge"

# =================================================================================================
#     Sample Names
# =================================================================================================

# The pileup format does not contain sample names. So instead, we produce list files for
# the mpileup files that contain multiple columns of data, where the list file contains
# the sample names, one per line. This can then for example also be used as input for the
# `--sample-name-list` input option in grenedalf!

# Attention: This order needs to be the same as what the below pilup rules use,
# otherwise we miss the point of writing these lists!

rule mpileup_all_sample_names:
    output:
        "mpileup/all-sample-names.txt"
    run:
        with open(output[0], "w") as out:
            for sample in config["global"]["sample-names"]:
                out.write(sample + "\n")

localrules: mpileup_all_sample_names

# =================================================================================================
#     Pileup
# =================================================================================================

# We offer several ways to create pileup files, depending on what the user configured:
# Either all in one file, all merged, or each sample as an individual file.
# We always use the merged units per sample here. This is a change after grenepipe v0.11.1,
# as it did not seem like treating each unit of each sample as an individual file made much sense.

rule mpileup_individual_sample:
    input:
        bam=get_sample_bams_wildcards, # provided in mapping.smk
        reference_genome=config["data"]["reference-genome"]
    output:
        "mpileup/{sample}.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/samples-{sample}.log"
    wrapper:
        "v3.13.6/bio/samtools/mpileup"

rule mpileup_all_samples:
    input:
        bam=get_all_bams(),
        reference_genome=config["data"]["reference-genome"],
        list="mpileup/all-sample-names.txt"
    output:
        "mpileup/all-samples.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/all-merged-units.log"
    wrapper:
        "v3.13.6/bio/samtools/mpileup"

rule mpileup_all_merged_samples:
    input:
        bam="mpileup/all.merged.bam",
        reference_genome=config["data"]["reference-genome"]
    output:
        "mpileup/all-merged-samples.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/all-merged-samples.log"
    wrapper:
        "v3.13.6/bio/samtools/mpileup"

# =================================================================================================
#     All pileups, but not SNP calling
# =================================================================================================

# This alternative target rule executes all steps up to th mapping and mpileup creation.
# This is the same as the above all_bams rule, but additionally also requests the pileups.
rule all_pileups:
    input:
        expand(
            "mpileup/{sample}.mpileup.gz",
            sample=config["global"]["sample-names"]
        ) if "individual-samples" in config["settings"]["pileups"] else [],
        (
            "mpileup/all-samples.mpileup.gz"
            if "all-samples" in config["settings"]["pileups"]
            else []
        ),
        (
            "mpileup/all-merged-samples.mpileup.gz"
            if "all-merged-samples" in config["settings"]["pileups"]
            else []
        )
    output:
        touch("mpileup/all-pileups.done")

# The `all_pileups` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all_pileups
