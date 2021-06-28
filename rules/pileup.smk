# =================================================================================================
#     Merge bams
# =================================================================================================

# Simple wildcard resolution.
def get_sample_bams_wildcards(wildcards):
    return get_sample_bams(wildcards.sample)

rule samtools_merge_units:
    input:
        get_sample_bams_wildcards
    output:
        temp("mpileup/{sample}.merged.bam")
    params:
        # Need to set the file overwrite flag for pipes to work here,
        # see https://github.com/samtools/samtools/issues/1437
        config["params"]["samtools"]["merge"] + " -f"
    threads:
        # Samtools takes additional threads through its option -@
        # This value - 1 will be sent to -@
        config["params"]["samtools"]["merge-threads"]
    log:
        "logs/samtools/mpileup/merge-{sample}.log"
    wrapper:
        "0.74.0/bio/samtools/merge"

rule samtools_merge_all:
    input:
        get_all_bams()
    output:
        pipe("mpileup/all.merged.bam")
    params:
        # Need file overwrite flag, see above.
        config["params"]["samtools"]["merge"] + " -f"
    threads:
        # Samtools takes additional threads through its option -@
        # This value - 1 will be sent to -@
        config["params"]["samtools"]["merge-threads"]
    log:
        "logs/samtools/mpileup/merge-all.log"
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

rule mpileup_samples_individual_units:
    input:
        bam=get_sample_bams_wildcards,
        reference_genome=config["data"]["reference"]["genome"]
    output:
        "mpileup/{sample}-individual-units.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/samples-individual-units-{sample}.log"
    wrapper:
        "0.74.0/bio/samtools/mpileup"

rule mpileup_samples_merged_units:
    input:
        bam="mpileup/{sample}.merged.bam",
        reference_genome=config["data"]["reference"]["genome"]
    output:
        "mpileup/{sample}-merged-units.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/samples-merged-units-{sample}.log"
    wrapper:
        "0.74.0/bio/samtools/mpileup"

rule mpileup_all_individual_units:
    input:
        bam=get_all_bams(),
        reference_genome=config["data"]["reference"]["genome"]
    output:
        "mpileup/all-individual-units.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/all-individual-units.log"
    wrapper:
        "0.74.0/bio/samtools/mpileup"

rule mpileup_all_merged_units:
    input:
        bam=expand(
            "mpileup/{sample}.merged.bam",
            sample=config["global"]["sample-names"]
        ),
        reference_genome=config["data"]["reference"]["genome"]
    output:
        "mpileup/all-merged-units.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/all-merged-units.log"
    wrapper:
        "0.74.0/bio/samtools/mpileup"

rule mpileup_all_merged_samples:
    input:
        bam="mpileup/all.merged.bam",
        reference_genome=config["data"]["reference"]["genome"]
    output:
        "mpileup/all-merged-samples.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/all-merged-samples.log"
    wrapper:
        "0.74.0/bio/samtools/mpileup"
