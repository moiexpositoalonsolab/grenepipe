# =================================================================================================
#     Merge bams
# =================================================================================================

# Simple wildcard resolution.
def get_sample_bams_wildcards(wildcards):
    return get_sample_bams(wildcards.sample)

rule samtools_merge_all:
    input:
        get_all_bams()
    output:
        pipe("mpileup/all.merged.bam")
    params:
        # Need to set the file overwrite flag for pipes to work here,
        # see https://github.com/samtools/samtools/issues/1437
        config["params"]["samtools"]["merge"] + " -f"
    threads:
        # Samtools takes additional threads through its option -@
        # This value - 1 will be sent to -@
        config["params"]["samtools"]["merge-threads"]
    wrapper:
        "0.74.0/bio/samtools/merge"

rule samtools_merge_units:
    input:
        get_sample_bams_wildcards
    output:
        pipe("mpileup/{sample}.merged.bam")
    params:
        # Need file overwrite flag, see above.
        config["params"]["samtools"]["merge"] + " -f"
    threads:
        # Samtools takes additional threads through its option -@
        # This value - 1 will be sent to -@
        config["params"]["samtools"]["merge-threads"]
    wrapper:
        "0.74.0/bio/samtools/merge"

# =================================================================================================
#     Pileup
# =================================================================================================

# We offer several ways to create pileup files, depending on what the user configured:
# Either all in one file, or each sample as an individual file, and for both these choices,
# either with everything merged, or as individual parts of the resulting file(s).
# That is, all-merged would be a single file with everything piled up, all-individual
# would be an mpileup with all samples and all units separately

def get_mpileup_all_bams():
    if config["settings"]["pileup-merge"] == "none":
        return get_all_bams()
    elif config["settings"]["pileup-merge"] == "units":
        return expand(
            "mpileup/{sample}.merged.bam",
            sample=config["global"]["sample-names"]
        )
    elif config["settings"]["pileup-merge"] == "all":
        return "mpileup/all.merged.bam"
    else:
        raise Exception(
            "Invalid value for settings:pileup-merge: " +  config["settings"]["pileup-merge"]
        )

rule mpileup_all:
    input:
        bam=get_mpileup_all_bams(),
        reference_genome=config["data"]["reference"]["genome"]
    output:
        "mpileup/all.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/all.log"
    wrapper:
        "0.74.0/bio/samtools/mpileup"

def get_mpileup_sample_bams(wildcards):
    if config["settings"]["pileup-merge"] == "none":
        return get_sample_bams_wildcards(wildcards)
    elif config["settings"]["pileup-merge"] == "units":
        return "mpileup/{}.merged.bam".format( wildcards.sample )
    elif config["settings"]["pileup-merge"] == "all":
        raise Exception(
            "Invalid value for settings:pileup-merge: all. "
            "Cannot use this value when creating per-sample mpileup files."
        )
    else:
        raise Exception(
            "Invalid value for settings:pileup-merge: " +  config["settings"]["pileup-merge"]
        )

rule mpileup_sample:
    input:
        bam=get_mpileup_sample_bams,
        reference_genome=config["data"]["reference"]["genome"]
    output:
        "mpileup/{sample}.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/{sample}.log"
    wrapper:
        "0.74.0/bio/samtools/mpileup"
