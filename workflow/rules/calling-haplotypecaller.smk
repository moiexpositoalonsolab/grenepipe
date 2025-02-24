import platform

# =================================================================================================
#     Variant Calling
# =================================================================================================


# Depending on the config, we need to either specify the restricted regions, or the combined
# contigs as the intervals in which to call. If neither is given, we call simply per contig
# in the ref genome.
def get_gatk_intervals(wildcards):
    if config["settings"].get("restrict-regions"):
        return "calling/regions/{}.bed".format(wildcards.contig)
    if config["settings"].get("contig-group-size"):
        return "calling/contig-groups/{}.bed".format(wildcards.contig)
    return wildcards.contig


# We also need a function that returns the file or an empty list, so that we can use this
# in the inout of rules, in order to ensure that the files are present.
def get_gatk_interval_files(wildcards):
    if config["settings"].get("restrict-regions"):
        return "calling/regions/{}.bed".format(wildcards.contig)
    if config["settings"].get("contig-group-size"):
        return "calling/contig-groups/{}.bed".format(wildcards.contig)
    return []


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
    log:
        "logs/calling/gatk-haplotypecaller/{sample}-{contig}.log",
    benchmark:
        "benchmarks/calling/gatk-haplotypecaller/{sample}-{contig}.log"
    # Need to set threads here so that snakemake can plan the job scheduling properly
    threads: config["params"]["gatk"]["HaplotypeCaller-threads"]
    params:
        # The intervals param here is where the contig variable is propagated to haplotypecaller.
        # Contigs are used as long as no restrict-regions are given in the config file.
        intervals=get_gatk_intervals,
        extra=config["params"]["gatk"].get("HaplotypeCaller-extra", ""),
        java_opts=config["params"]["gatk"].get("HaplotypeCaller-java-opts", ""),
    resources:
        # Increase time limit in factors of 24h, if the job fails due to time limit.
        # time = lambda wildcards, input, threads, attempt: int(1440 * int(attempt))
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


# Recommended way of GATK to combine GVCFs these days.
rule genomics_db_import:
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
            "calling/called/{sample}-{{contig}}.g.vcf.gz", sample=config["global"]["sample-names"]
        ),
        indices=expand(
            "calling/called/{sample}-{{contig}}.g.vcf.gz.tbi",
            sample=config["global"]["sample-names"],
        ),
        done=expand(
            "calling/called/{sample}-{{contig}}.g.vcf.gz.done",
            sample=config["global"]["sample-names"],
        ),
        # Same as above, we need a dummy for the intervals to ensure the files are present.
        intervals_dummy=get_gatk_interval_files,
    output:
        db=directory("calling/genomics_db/{contig}"),
        done=touch("calling/genomics_db/{contig}.done"),
    log:
        "logs/calling/gatk-genomicsdbimport/{contig}.log",
    benchmark:
        "benchmarks/calling/gatk-genomicsdbimport/{contig}.log"
    params:
        # Here, we actually use the intervals to provide them to the wrapper.
        intervals=get_gatk_intervals,
        db_action="create",
        extra=" --reference "
        + config["data"]["reference-genome"]
        + " --sequence-dictionary "
        + genome_dict()
        + " "
        + config["params"]["gatk"].get("GenomicsDBImport-extra", ""),
        java_opts=config["params"]["gatk"].get("GenomicsDBImport-java-opts", ""),
    # threads: 2
    resources:
        mem_mb=config["params"]["gatk"].get("GenomicsDBImport-mem-mb", 1024),
        tmpdir=config["params"]["gatk"].get("GenomicsDBImport-temp-dir", ""),
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/genomicsdbimport"


# Old way, using GATK CombineGVCFs, which is slow when run on many samples.
# We still offer it for compatibility and completeness, but recommend using GenomicsDBImport instead.
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
            "calling/called/{sample}-{{contig}}.g.vcf.gz", sample=config["global"]["sample-names"]
        ),
        indices=expand(
            "calling/called/{sample}-{{contig}}.g.vcf.gz.tbi",
            sample=config["global"]["sample-names"],
        ),
        done=expand(
            "calling/called/{sample}-{{contig}}.g.vcf.gz.done",
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
        "benchmarks/calling/gatk-combine-gvcfs/{contig}.log"
    # group:
    #     "gatk_calls_combine"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/combinegvcfs"


# =================================================================================================
#     Genotype Variants
# =================================================================================================


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
        # Get the GVCF or GenomicsDB input, depending on which tool is requested in the config.
        gvcf=(
            "calling/combined/all.{contig}.g.vcf.gz"
            if not config["params"]["gatk"].get("use-GenomicsDBImport", True)
            else []
        ),
        gvcf_done=(
            "calling/combined/all.{contig}.g.vcf.gz.done"
            if not config["params"]["gatk"].get("use-GenomicsDBImport", True)
            else []
        ),
        genomicsdb=(
            "calling/genomics_db/{contig}"
            if config["params"]["gatk"].get("use-GenomicsDBImport", True)
            else []
        ),
        genomicsdb_done=(
            "calling/genomics_db/{contig}.done"
            if config["params"]["gatk"].get("use-GenomicsDBImport", True)
            else []
        ),
        # If known variants are set in the config, use them, and require the index file as well.
        known=config["data"]["known-variants"],
        knownidx=(
            config["data"]["known-variants"] + ".tbi" if config["data"]["known-variants"] else []
        ),
        # Same as above, we need a dummy for the intervals to ensure the files are present.
        intervals_dummy=get_gatk_interval_files,
    output:
        vcf=(
            "calling/genotyped/all.{contig}.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/genotyped/all.{contig}.vcf.gz")
        ),
        done=touch("calling/genotyped/all.{contig}.vcf.gz.done"),
    params:
        # Again, we here use the intervals to provide them to the wrapper.
        intervals=get_gatk_intervals,
        extra=" --sequence-dictionary "
        + genome_dict()
        + " "
        + config["params"]["gatk"]["GenotypeGVCFs-extra"],
        java_opts=config["params"]["gatk"]["GenotypeGVCFs-java-opts"],
    resources:
        mem_mb=config["params"]["gatk"].get("GenotypeGVCFs-mem-mb", 1024),
    log:
        "logs/calling/gatk-genotype-gvcfs/{contig}.log",
    benchmark:
        "benchmarks/calling/gatk-genotype-gvcfs/{contig}.log"
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
        java_opts=config["params"]["picard"].get("MergeVcfs-java-opts", ""),
        extra=(
            " --USE_JDK_DEFLATER true --USE_JDK_INFLATER true"
            if platform.system() == "Darwin"
            else ""
        ),
    resources:
        mem_mb=config["params"]["picard"].get("MergeVcfs-mem-mb", 1024),
    log:
        "logs/calling/picard-merge-genotyped.log",
    benchmark:
        "benchmarks/calling/picard-merge-genotyped.log"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "v5.7.0/bio/picard/mergevcfs"
