# =================================================================================================
#     Variant Selection Helper
# =================================================================================================

def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format( "SNP" if wildcards.vartype == "snvs" else "INDEL" )

rule select_calls:
    input:
        ref=config["data"]["reference"]["genome"],
        vcf=config["rundir"] + "genotyped/all.vcf.gz",

        # bcftools does not automatically create vcf index files, so we need to specifically request them...
        # ... but the picard merge tool that we use right now does create tbi files, so all good atm.
        # tbi=config["rundir"] + "genotyped/all.vcf.gz.tbi" if config["settings"]["calling-tool"] == "bcftools" else []
    output:
        vcf=temp(config["rundir"] + "filtered/all.select-{vartype}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        config["rundir"] + "logs/gatk/selectvariants/{vartype}.log"
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
        vcf=config["rundir"] + "filtered/all.select-{vartype}.vcf.gz"
    output:
        vcf=temp(config["rundir"] + "filtered/all.{vartype}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    log:
        config["rundir"] + "logs/gatk/variantfiltration/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/variantfiltration"

# Machine learning based recalibration of quality scores instead of hard filtering.
# To use this, set the settings:vqsr switch to true in the config.yaml
rule recalibrate_calls:
    input:
        vcf=config["rundir"] + "filtered/all.select-{vartype}.vcf.gz"
    output:
        vcf=temp(config["rundir"] + "filtered/all.{vartype}.recalibrated.vcf.gz")
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]
    log:
        config["rundir"] + "logs/gatk/variantrecalibrator/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/variantrecalibrator"

# =================================================================================================
#     Merge Filtered Variants
# =================================================================================================

rule merge_calls:
    input:
        vcf=expand(config["rundir"] + "filtered/all.{vartype}.{filtertype}.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated" if config["settings"]["vqsr"] else "hardfiltered")
    output:
        vcf=protected(config["rundir"] + "filtered/all.vcf.gz")
    log:
        config["rundir"] + "logs/picard/merge-filtered.log"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"
