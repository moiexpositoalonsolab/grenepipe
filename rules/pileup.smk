# =================================================================================================
#     Merge bams
# =================================================================================================

# This per-sample merging is also present in the HAF-pipe rules in frequency.smk,
# and in the qualimap rule in qc.smk; meaning that this step is potentially executed
# multiple times independently, for different purposes. For now, we keep it this way,
# to have more control over the specifics, but we might change that in the future.

rule mpileup_merge_unit_bams:
    input:
        get_sample_bams_wildcards # provided in mapping.smk
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

rule mpileup_merge_all:
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
#     Sample Names
# =================================================================================================

# The pileup format does not contain sample names. So instead, we produce list files for
# the mpileup files that contain multiple columns of data, where the list file contains
# the sample names, one per line. This can then for example also be used as input for the
# `--sample-name-list` input option in grenedalf!

# Attention: This order needs to be the same as what the below pilup rules use,
# otherwise we miss the point of writing these lists!

rule mpileup_all_individual_units_names:
    output:
        "mpileup/all-individual-units.names.txt"
    run:
        with open(output[0], "w") as out:
            for su in config["global"]["sample-units"]:
                out.write(su[0] + "-" + su[1] + "\n")

rule mpileup_all_merged_units_names:
    output:
        "mpileup/all-merged-units.names.txt"
    run:
        with open(output[0], "w") as out:
            for sample in config["global"]["sample-names"]:
                out.write(sample + "\n")

localrules: mpileup_all_individual_units_names
localrules: mpileup_all_merged_units_names

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
        bam=get_sample_bams_wildcards, # provided in mapping.smk
        reference_genome=config["data"]["reference-genome"]
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
        reference_genome=config["data"]["reference-genome"]
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
        reference_genome=config["data"]["reference-genome"],
        list="mpileup/all-individual-units.names.txt"
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
        reference_genome=config["data"]["reference-genome"],
        list="mpileup/all-merged-units.names.txt"
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
        reference_genome=config["data"]["reference-genome"]
    output:
        "mpileup/all-merged-samples.mpileup.gz"
    params:
        extra=config["params"]["samtools"]["pileup"]
    log:
        "logs/samtools/mpileup/all-merged-samples.log"
    wrapper:
        "0.74.0/bio/samtools/mpileup"

# =================================================================================================
#     All pileups, but not SNP calling
# =================================================================================================

# This alternative target rule executes all steps up to th mapping and mpileup creation.
# This is the same as the above all_bams rule, but additionally also requests the pileups.
rule all_pileups:
    input:
        expand(
            "mpileup/{sample}-individual-units.mpileup.gz",
            sample=config["global"]["sample-names"]
        ) if "samples-individual-units" in config["settings"]["pileups"] else [],
        expand(
            "mpileup/{sample}-merged-units.mpileup.gz",
            sample=config["global"]["sample-names"]
        ) if "samples-merged-units" in config["settings"]["pileups"] else [],
        "mpileup/all-individual-units.mpileup.gz" if "all-individual-units" in config["settings"]["pileups"] else [],
        "mpileup/all-merged-units.mpileup.gz" if "all-merged-units" in config["settings"]["pileups"] else [],
        "mpileup/all-merged-samples.mpileup.gz" if "all-merged-samples" in config["settings"]["pileups"] else [],

# The `all_pileups` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all_pileups
