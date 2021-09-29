# =================================================================================================
#     Variant Selection Helper
# =================================================================================================

def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format( "SNP" if wildcards.vartype == "snvs" else "INDEL" )

rule select_calls:
    input:
        ref=config["data"]["reference"]["genome"],
        vcf="genotyped/all.vcf.gz",
        refdict=genome_dict(),

        # bcftools does not automatically create vcf index files, so we need to specifically request them...
        # ... but the picard merge tool that we use right now does create tbi files, so all good atm.
        # tbi="genotyped/all.vcf.gz.tbi" if config["settings"]["calling-tool"] == "bcftools" else []
    output:
        vcf=temp("filtered/all.select-{vartype}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        "logs/gatk/selectvariants/{vartype}.log"
    benchmark:
        "benchmarks/gatk/selectvariants/{vartype}.bench.log"
    group:
        "filtering"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "0.27.1/bio/gatk/selectvariants"

# =================================================================================================
#     Filtering
# =================================================================================================

def get_filter(wildcards):
    return { "snv-hard-filter": config["params"]["variantfiltration-hard"][wildcards.vartype] }

# Simple hard filter, used by default
rule hard_filter_calls:
    input:
        ref=config["data"]["reference"]["genome"],
        vcf="filtered/all.select-{vartype}.vcf.gz"
    output:
        vcf=temp("filtered/all.{vartype}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    log:
        "logs/gatk/variantfiltration/{vartype}.log"
    benchmark:
        "benchmarks/gatk/variantfiltration/{vartype}.bench.log"
    group:
        "filtering"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "0.27.1/bio/gatk/variantfiltration"

# Machine learning based recalibration of quality scores instead of hard filtering.
# To use this, there was a config settings:vqsr switch, but re removed this now,
# as this method is not supported at the moment. We keep the code here for reference,
# but it never worked properly. If it becomes necessary in the future to use this method,
# we probably need to re-implement this rule completely.
# rule recalibrate_calls:
#     input:
#         vcf="filtered/all.select-{vartype}.vcf.gz"
#     output:
#         vcf=temp("filtered/all.{vartype}.recalibrated.vcf.gz")
#     params:
#         extra=config["params"]["gatk"]["VariantRecalibrator"]
#     log:
#         "logs/gatk/variantrecalibrator/{vartype}.log"
#     benchmark:
#         "benchmarks/gatk/variantrecalibrator/{vartype}.bench.log"
#     group:
#         "filtering"
#     wrapper:
#         "0.27.1/bio/gatk/variantrecalibrator"

# =================================================================================================
#     Merge Filtered Variants
# =================================================================================================

rule merge_calls:
    input:
        vcf=expand("filtered/all.{vartype}.{filtertype}.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated" if False else "hardfiltered")
                   # Keeping vqsr here for reference. Not used at the moment.
                   # filtertype="recalibrated" if config["settings"]["vqsr"] else "hardfiltered")
    output:
        vcf=protected("filtered/all.vcf.gz")
    log:
        "logs/picard/merge-filtered.log"
    benchmark:
        "benchmarks/picard/merge-filtered.bench.log"
    group:
        "filtering"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"
