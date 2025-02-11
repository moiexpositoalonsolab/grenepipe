import platform

# =================================================================================================
#     Variant Calling
# =================================================================================================


# Combine all params to call gatk. We may want to set regions, we set that bit of multithreading
# that gatk is capable of (not much, but the best we can do without spark...), and we add
# all additional params from the config file.
def get_gatk_call_variants_params(wildcards, input):
    return (
        get_gatk_regions_param(
            regions=input.regions, default="--intervals '{}'".format(wildcards.contig)
        )
        + " --native-pair-hmm-threads "
        + str(config["params"]["gatk"]["HaplotypeCaller-threads"])
        + " "
        + config["params"]["gatk"]["HaplotypeCaller-extra"]
    )


rule call_variants:
    input:
        # Get the sample data.
        bam=get_sample_bams_wildcards,
        bai=get_sample_bais_wildcards,
        done=get_sample_bams_wildcards_done,
        # Get the reference genome, as well as its indices.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        refdict=genome_dict(),
        # If known variants are set in the config, use then, and require the index file as well.
        known=config["data"]["known-variants"],
        knownidx=(
            config["data"]["known-variants"] + ".tbi" if config["data"]["known-variants"] else []
        ),
        # Further settings for region constraint filter.
        # regions="calling/regions/{contig}.bed" if config["settings"].get("restrict-regions") else []
        regions=(
            "calling/regions/{contig}.bed"
            if (config["settings"].get("restrict-regions"))
            else (
                "calling/contig-groups/{contig}.bed"
                if (config["settings"].get("contig-group-size"))
                else []
            )
        ),
    output:
        gvcf=(
            "calling/called/{sample}.{contig}.g.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/called/{sample}.{contig}.g.vcf.gz")
        ),
        # gvcf=protected("calling/called/{sample}.{contig}.g.vcf.gz")
        gtbi=(
            "calling/called/{sample}.{contig}.g.vcf.gz.tbi"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/called/{sample}.{contig}.g.vcf.gz.tbi")
        ),
        done=touch("calling/called/{sample}.{contig}.g.vcf.gz.done"),
    log:
        "logs/calling/gatk-haplotypecaller/{sample}.{contig}.log",
    benchmark:
        "benchmarks/calling/called/gatk-haplotypecaller/{sample}.{contig}.log"
    # Need to set threads here so that snakemake can plan the job scheduling properly
    threads: config["params"]["gatk"]["HaplotypeCaller-threads"]
    # resources:
    # Increase time limit in factors of 24h, if the job fails due to time limit.
    # time = lambda wildcards, input, threads, attempt: int(1440 * int(attempt))
    params:
        # The function here is where the contig variable is propagated to haplotypecaller.
        # Took me a while to figure this one out...
        # Contigs are used as long as no restrict-regions are given in the config file.
        extra=get_gatk_call_variants_params,
        java_opts=config["params"]["gatk"]["HaplotypeCaller-java-opts"],
    resources:
        mem_mb=config["params"]["gatk"].get("HaplotypeCaller-mem-mb", 1024),
    group:
        "call_variants"
    conda:
        # Need to specify, yet again...
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/haplotypecaller"


# Deactivated the below, as this was causing trouble. Got the warning
#     Warning: the following output files of rule vcf_index_gatk were not present when the DAG was created:
#     {'called/S3.chloroplast.g.vcf.gz.tbi'}
# for all files, indicating that the above rule indeed does produce them.
# However, having an extra rule for that caused that rule to _sometimes_ be executed, so that
# the tbi file would have a later time stamp, and it seems likely that this then caused other
# rules to want to update as well, meaning that the snp calling was repeated?!
# I hope that this fix this problem...

# # Stupid GATK sometimes writes out index files, and sometimes not, and it is not clear at all
# # when that is happening and when not. Let's try with a rule, and see if it works even if the file
# # is present sometimes... hopefully snakemake is smart enough for that.
# rule vcf_index_gatk:
#     input:
#         "calling/{file}.g.vcf.gz"
#     output:
#         "calling/{file}.g.vcf.gz.tbi"
#     params:
#         # pass arguments to tabix (e.g. index a vcf)
#         "-p vcf"
#     log:
#         "logs/tabix/{file}.log"
#     group:
#         "call_variants"
#     wrapper:
#         "0.55.1/bio/tabix"

# =================================================================================================
#     Combining Calls
# =================================================================================================


rule combine_calls:
    input:
        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        refdict=genome_dict(),
        # Get the sample data, including indices.
        gvcfs=expand(
            "calling/called/{sample}.{{contig}}.g.vcf.gz", sample=config["global"]["sample-names"]
        ),
        indices=expand(
            "calling/called/{sample}.{{contig}}.g.vcf.gz.tbi",
            sample=config["global"]["sample-names"],
        ),
        done=expand(
            "calling/called/{sample}.{{contig}}.g.vcf.gz.done",
            sample=config["global"]["sample-names"],
        ),
    output:
        gvcf=(
            "calling/combined/all.{contig}.g.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/combined/all.{contig}.g.vcf.gz")
        ),
        done=touch("calling/combined/all.{contig}.g.vcf.gz.done"),
    params:
        extra=config["params"]["gatk"]["CombineGVCFs-extra"]
        + (
            " --dbsnp " + config["data"]["known-variants"] + " "
            if config["data"]["known-variants"]
            else ""
        ),
        java_opts=config["params"]["gatk"]["CombineGVCFs-java-opts"],
    resources:
        mem_mb=config["params"]["gatk"].get("CombineGVCFs-mem-mb", 1024),
    log:
        "logs/calling/gatk-combine-gvcfs/{contig}.log",
    benchmark:
        "benchmarks/calling/called/gatk-combine-gvcfs/{contig}.log"
    # group:
    #     "gatk_calls_combine"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        refdict=genome_dict(),
        gvcf="calling/combined/all.{contig}.g.vcf.gz",
        done="calling/combined/all.{contig}.g.vcf.gz.done",
    output:
        vcf=(
            "calling/genotyped/all.{contig}.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/genotyped/all.{contig}.vcf.gz")
        ),
        done=touch("calling/genotyped/all.{contig}.vcf.gz.done"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs-extra"]
        + (
            " --dbsnp " + config["data"]["known-variants"] + " "
            if config["data"]["known-variants"]
            else ""
        ),
        java_opts=config["params"]["gatk"]["GenotypeGVCFs-java-opts"],
    resources:
        mem_mb=config["params"]["gatk"].get("GenotypeGVCFs-mem-mb", 1024),
    log:
        "logs/calling/gatk-genotype-gvcfs/{contig}.log",
    benchmark:
        "benchmarks/calling/genotyped/gatk-genotype-gvcfs/{contig}.log"
    # group:
    #     "gatk_calls_combine"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/genotypegvcfs"


# =================================================================================================
#     Merging Variants
# =================================================================================================


# Need an input function to work with the fai checkpoint
def merge_variants_vcfs_input(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("calling/genotyped/all.{contig}.vcf.gz", contig=get_contigs(fai))


# Also need to trick snakemake into completing all files...
def merge_variants_vcfs_input_done(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("calling/genotyped/all.{contig}.vcf.gz.done", contig=get_contigs(fai))


rule merge_variants:
    input:
        # fai is needed to calculate aggregation over contigs below.
        # This is the step where the genome is split into its contigs for parallel execution.
        # The get_fai() function uses a snakemake checkpoint to make sure that the fai is
        # produced before we use it here to get its content.
        ref=get_fai,
        contig_groups=contigs_groups_input,
        # vcfs=lambda w: expand("calling/genotyped/all.{contig}.vcf.gz", contig=get_contigs())
        vcfs=merge_variants_vcfs_input,
        done=merge_variants_vcfs_input_done,
    output:
        vcf="calling/genotyped-all.vcf.gz",
        done=touch("calling/genotyped-all.vcf.gz.done"),
    params:
        # See duplicates-picard.smk for the reason whe need this on MacOS.
        java_opts=config["params"]["picard"]["MergeVcfs-java-opts"],
        extra=(
            " --USE_JDK_DEFLATER true --USE_JDK_INFLATER true" if platform.system() == "Darwin" else ""
        ),
    resources:
        mem_mb=config["params"]["picard"].get("MergeVcfs-mem-mb", 1024),
    log:
        "logs/calling/picard-merge-genotyped.log",
    benchmark:
        "benchmarks/calling/genotyped/picard/merge-genotyped.log"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "v5.7.0/bio/picard/mergevcfs"
