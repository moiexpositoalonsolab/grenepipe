import platform

# =================================================================================================
#     Helpers
# =================================================================================================


# Depending on the config, we need to either specify the restricted regions, or the combined
# contigs as the intervals in which to call. If neither is given, we call simply per contig
# in the ref genome.
def get_gatk_intervals(wildcards):
    if config["settings"].get("restrict-regions"):
        return "calling/regions/{}.bed".format(wildcards.contig)
    if config["settings"].get("contig-group-size", 0):
        return "calling/contig-groups/{}.bed".format(wildcards.contig)
    return wildcards.contig


# We also need a function that returns the file or an empty list, so that we can use this
# in the inout of rules, in order to ensure that the files are present.
def get_gatk_interval_files(wildcards):
    if config["settings"].get("restrict-regions"):
        return "calling/regions/{}.bed".format(wildcards.contig)
    if config["settings"].get("contig-group-size", 0):
        return "calling/contig-groups/{}.bed".format(wildcards.contig)
    return []


# =================================================================================================
#     Haplotype Calling
# =================================================================================================


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
        # We need this here as an unused input, so that the bed files are guaranteed to be created
        # beforehand if needed. They are however actually provided to the wrapper via the params.
        # We cannot provide them here, in the default case, it's a contig name, which is not a file.
        # Hence, this is only relevant with restrict regions or contig groups.
        intervals_dummy=get_gatk_interval_files,
    output:
        gvcf=(
            "calling/called/{sample}-{contig}.g.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/called/{sample}-{contig}.g.vcf.gz")
        ),
        # gvcf=protected("calling/called/{sample}-{contig}.g.vcf.gz")
        gtbi=(
            "calling/called/{sample}-{contig}.g.vcf.gz.tbi"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/called/{sample}-{contig}.g.vcf.gz.tbi")
        ),
        done=touch("calling/called/{sample}-{contig}.g.vcf.gz.done"),
    params:
        # The intervals param here is where the contig variable is propagated to haplotypecaller.
        # Contigs are used as long as no restrict-regions are given in the config file.
        intervals=get_gatk_intervals,
        extra=config["params"]["gatk"].get("HaplotypeCaller-extra", ""),
        java_opts=config["params"]["gatk"].get("HaplotypeCaller-java-opts", ""),
    threads: get_rule_threads("call_variants")
    log:
        "logs/calling/gatk-haplotypecaller/{sample}-{contig}.log",
    benchmark:
        "benchmarks/calling/gatk-haplotypecaller/{sample}-{contig}.log"
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
#     Genotype Calling
# =================================================================================================

# Combining the per-sample data into a format that the joint
# variant calling can work from, using different stragegies.
# For each, we have a separate implementation file, to keep things managed.
# They start with the above haplotype calles, and end with per-contig vcfs
# that are then merged below.
if not config["params"]["gatk"].get("use-GenomicsDBImport", True):

    include: "calling-gatk-genotype-gvcfs.smk"

elif config["params"]["gatk"].get("GenomicsDBImport-interval-size", 0) == 0:

    include: "calling-gatk-genotype-db-all.smk"

else:

    include: "calling-gatk-genotype-db-split.smk"


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
        java_opts=config["params"]["picard"].get("MergeVcfs-java-opts", ""),
        extra=(
            " --USE_JDK_DEFLATER true --USE_JDK_INFLATER true"
            if platform.system() == "Darwin"
            else ""
        ),
    threads: get_rule_threads("merge_variants")
    log:
        "logs/calling/picard-merge-genotyped.log",
    benchmark:
        "benchmarks/calling/picard-merge-genotyped.log"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "v5.7.0/bio/picard/mergevcfs"
