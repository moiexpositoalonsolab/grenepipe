# =================================================================================================
#     Restrict Regions
# =================================================================================================

if "restrict-regions" in config["data"]:
    rule compose_regions:
        input:
            config["data"]["reference"].get("restrict-regions")
        output:
            "called/{contig}.regions.bed"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"

# =================================================================================================
#     Variant Calling
# =================================================================================================

def get_sample_bams(wildcards):
    """
    Get all aligned reads of given sample, with all its units.
    This is where all units are merged together. The function also automatically gets
    which of the mapping resutls to use, depending on the config setting (whether to remove
    duplicates, and whether to recalibrate the base qualities), by using the get_mapping_result
    function, that gives the respective files depending on the config.
    """
    return expand(get_mapping_result(), sample=wildcards.sample, unit=samples.loc[wildcards.sample].unit)

rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["data"]["reference"]["genome"],
        known=config["data"]["reference"].get("known-variants"), # empty if key not present
        regions="called/{contig}.regions.bed" if config["data"]["reference"].get("restrict-regions") else []
    output:
        gvcf=protected("called/{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log"
    params:
        # The function here is where the contig variable is propagated to haplotypecaller.
        # Took me a while to figure this one out...
        extra=get_call_variants_params
    wrapper:
        "0.51.3/bio/gatk/haplotypecaller"

# =================================================================================================
#     Combining Calls
# =================================================================================================

rule combine_calls:
    input:
        ref=config["data"]["reference"]["genome"],
        gvcfs=expand("called/{sample}.{{contig}}.g.vcf.gz", sample=sample_indices)
    output:
        gvcf="called/all.{contig}.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.{contig}.log"
    wrapper:
        "0.51.3/bio/gatk/combinegvcfs"

rule genotype_variants:
    input:
        ref=config["data"]["reference"]["genome"],
        gvcf="called/all.{contig}.g.vcf.gz"
    output:
        vcf=temp("genotyped/all.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs.{contig}.log"
    wrapper:
        "0.51.3/bio/gatk/genotypegvcfs"

# =================================================================================================
#     Merging Variants
# =================================================================================================

def get_fai():
    return config["data"]["reference"]["genome"] + ".fai"

# contigs in reference genome
def get_contigs():
    return pd.read_table(get_fai(), header=None, usecols=[0], squeeze=True, dtype=str)

rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs=lambda w: expand("genotyped/all.{contig}.vcf.gz", contig=get_contigs()),
    output:
        vcf="genotyped/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    wrapper:
        "0.51.3/bio/picard/mergevcfs"
