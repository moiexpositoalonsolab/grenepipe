# =================================================================================================
#     Variant Calling
# =================================================================================================

# Return the bam file(s) for a given sample
def get_sample_bams_wildcard(wildcards):
    """
    Get all aligned reads of given sample, with all its units.
    This is where all units are merged together. The function also automatically gets
    which of the mapping resutls to use, depending on the config setting (whether to remove
    duplicates, and whether to recalibrate the base qualities), by using the get_mapping_result
    function, that gives the respective files depending on the config.
    """
    return expand(
        get_mapping_result(),
        sample=wildcards.sample,
        unit=config["global"]["samples"].loc[wildcards.sample].unit
    )

# Return the bai file(s) for a given sample
def get_sample_bais_wildcard(wildcards):
    return expand(
        get_mapping_result(True),
        sample=wildcards.sample,
        unit=config["global"]["samples"].loc[wildcards.sample].unit
    )

# Combine all params to call gatk. We may want to set regions, we set that bit of multithreading
# that gatk is capable of (not much, but the best we can do without spark...), and we add
# all additional params from the config file.
def get_gatk_call_variants_params(wildcards, input):
    return (
        get_gatk_regions_param(regions=input.regions, default="--intervals '{}'".format(wildcards.contig))
        + " --native-pair-hmm-threads " + str(config["params"]["gatk"]["HaplotypeCaller-threads"]) + " "
        + config["params"]["gatk"]["HaplotypeCaller"]
    )

rule call_variants:
    input:
        # Get the sample data.
        bam=get_sample_bams_wildcard,
        bai=get_sample_bais_wildcard,

        # Get the reference genome, as well as its indices.
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),
        refdict=genome_dict(),

        # If known variants are set in the config, use then, and require the index file as well.
        known=config["data"]["reference"].get("known-variants"), # empty if key not present
        knownidx=config["data"]["reference"]["known-variants"] + ".tbi" if config["data"]["reference"]["known-variants"] else [],

        # Further settings for region constraint filter.
        # regions="called/{contig}.regions.bed" if config["settings"].get("restrict-regions") else []
        regions="called/{contig}.regions.bed" if (
            config["settings"].get("restrict-regions")
        ) else (
            "contig-groups/{contig}.bed" if (
                config["settings"].get("contig-group-size")
            ) else []
        )
    output:
        gvcf=protected("called/{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log"
    benchmark:
        "benchmarks/gatk/haplotypecaller/{sample}.{contig}.bench.log"
    threads:
        # Need to set threads here so that snakemake can plan the job scheduling properly
        config["params"]["gatk"]["HaplotypeCaller-threads"]
    # resources:
        # Increase time limit in factors of 24h, if the job fails due to time limit.
        # time = lambda wildcards, input, threads, attempt: int(1440 * int(attempt))
    params:
        # The function here is where the contig variable is propagated to haplotypecaller.
        # Took me a while to figure this one out...
        # Contigs are used as long as no restrict-regions are given in the config file.
        extra=get_gatk_call_variants_params
    group:
        "call_variants"
    conda:
        # Need to specify, yet again...
        "../envs/gatk.yaml"
    shadow: "full"
    wrapper:
        "0.51.3/bio/gatk/haplotypecaller"

# Stupid GATK sometimes writes out index files, and sometimes not, and it is not clear at all
# when that is happening and when not. Let's try with a rule, and see if it works even if the file
# is present sometimes... hopefully snakemake is smart enough for that.
rule vcf_index_gatk:
    input:
        "called/{file}.g.vcf.gz"
    output:
        "called/{file}.g.vcf.gz.tbi"
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf"
    log:
        "logs/tabix/{file}.log"
    group:
        "call_variants"
    wrapper:
        "0.55.1/bio/tabix"

# =================================================================================================
#     Combining Calls
# =================================================================================================

rule combine_calls:
    input:
        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),

        # Get the sample data, including indices.
        gvcfs=expand(
            "called/{sample}.{{contig}}.g.vcf.gz",
            sample=config["global"]["sample-names"]
        ),
        indices=expand(
            "called/{sample}.{{contig}}.g.vcf.gz.tbi",
            sample=config["global"]["sample-names"]
        )
    output:
        gvcf=temp("called/all.{contig}.g.vcf.gz")
    log:
        "logs/gatk/combine-gvcfs/{contig}.log"
    benchmark:
        "benchmarks/gatk/combine-gvcfs/{contig}.bench.log"
    group:
        "gatk_calls_combine"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "0.51.3/bio/gatk/combinegvcfs"

rule genotype_variants:
    input:
        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),

        gvcf="called/all.{contig}.g.vcf.gz"
    output:
        vcf=temp("genotyped/all.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotype-gvcfs/{contig}.log"
    benchmark:
        "benchmarks/gatk/genotype-gvcfs/{contig}.bench.log"
    group:
        "gatk_calls_combine"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "0.51.3/bio/gatk/genotypegvcfs"

# =================================================================================================
#     Merging Variants
# =================================================================================================

# Need an input function to work with the fai checkpoint
def merge_variants_vcfs_input(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("genotyped/all.{contig}.vcf.gz", contig=get_contigs( fai ))

rule merge_variants:
    input:
        # fai is needed to calculate aggregation over contigs below.
        # This is the step where the genome is split into its contigs for parallel execution.
        # The get_fai() function uses a snakemake checkpoint to make sure that the fai is
        # produced before we use it here to get its content.
        ref=get_fai,

        # vcfs=lambda w: expand("genotyped/all.{contig}.vcf.gz", contig=get_contigs())
        vcfs=merge_variants_vcfs_input
    output:
        vcf="genotyped/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    benchmark:
        "benchmarks/picard/merge-genotyped.bench.log"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "0.51.3/bio/picard/mergevcfs"
