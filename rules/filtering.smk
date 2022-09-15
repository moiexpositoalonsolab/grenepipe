import json

# =================================================================================================
#     Variant Selection Helper
# =================================================================================================

rule select_calls:
    input:
        ref=config["data"]["reference-genome"],
        vcf="genotyped/all.vcf.gz",
        refdict=genome_dict(),

        # bcftools does not automatically create vcf index files, so we need to specifically request them...
        # ... but the picard merge tool that we use right now does create tbi files, so all good atm.
        # tbi="genotyped/all.vcf.gz.tbi" if config["settings"]["calling-tool"] == "bcftools" else []
    output:
        vcf=(
            "filtered/all.select-{vartype}.vcf.gz"
            if config["settings"]["keep-intermediate"]["filtering"]
            else temp("filtered/all.select-{vartype}.vcf.gz")
        )
    params:
        extra="--select-type-to-include {vartype}"
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
#     Hard Filtering
# =================================================================================================

def get_filter(wildcards):
    return {
        wildcards.vartype + "-hard-filter":
        config["params"]["variantfiltration-hard"][wildcards.vartype]
    }

# Simple hard filter, used by default
rule hard_filter_calls:
    input:
        ref=config["data"]["reference-genome"],
        vcf="filtered/all.select-{vartype}.vcf.gz"
    output:
        vcf=(
            "filtered/all.{vartype}.hardfiltered.vcf.gz"
            if config["settings"]["keep-intermediate"]["filtering"]
            else temp("filtered/all.{vartype}.hardfiltered.vcf.gz")
        )
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

# =================================================================================================
#     VQSR Recalibration and Filtering
# =================================================================================================

# We want our config file to be a bit flexible, and allow a dict or a list of dicts to be provided,
# so here we need to unpack those. Python wants this to be complicated, so we use a detour via json,
# see https://stackoverflow.com/a/27373027/4184258
def ordered_dict_to_dict(value):
    value = json.loads(json.dumps(value))
    if isinstance(value, list) and len(value) == 1 and isinstance(value[0], dict):
        value = value[0]
    if isinstance(value, dict):
        for k, v in value.items():
            if isinstance(v, list) and len(v) == 1 and isinstance(v[0], dict):
                v = v[0]
            if isinstance(v, dict):
                value[k] = ordered_dict_to_dict(v)
    # print(value)
    return value

# We need to turn our config file into named inputs for the resources,
# as this is the format that the below wrapper expects.
def get_variant_recalibrator_resource_files():
    rf = ordered_dict_to_dict( config["params"]["variantfiltration-vqsr"]["resource-files"] )
    if any(x in rf.keys() for x in ["ref", "vcf", "tbi"]):
        raise Exception("Cannot use resource named 'ref', 'vcf', or 'tbi' with our GATK VQSR.")
    return rf

# We also add the tbi files as dependencies, so that their existence is checked.
def get_variant_recalibrator_resource_file_tbis():
    rf = get_variant_recalibrator_resource_files()
    tbi = []
    for k in rf.keys():
        tbi.append( rf[k] + ".tbi" )
    return tbi

def get_variant_recalibrator_resources():
    res = ordered_dict_to_dict( config["params"]["variantfiltration-vqsr"]["resources"] )
    for rn, rv in res.items():
        if isinstance( rv, list ):
            nv = {}
            for erv in rv:
                for k,v in erv.items():
                    nv[ k ] = v
            res[rn] = nv
    return res

def get_variant_recalibrator_extra(wildcards):
    # Need to do this in a function, as wildcard replacement otherwise does not work within
    # config file dict retrieval. We additionally set the Rscript file, so that we get some plots.
    vartype = wildcards.vartype
    return (
        config["params"]["variantfiltration-vqsr"]["extra-variantrecalibrator-" + vartype] +
        " --rscript-file filtered/all." + vartype + ".vqsr-recal.plots.R"
    )

def get_apply_vqsr_extra(wildcards):
    # Same as above, need a function here for wildcard replacement
    return config["params"]["variantfiltration-vqsr"]["extra-applyvqsr-" + wildcards.vartype]

# Use the GATK VQSR machine learning based recalibration of quality scores instead of hard filtering.
# This is a two step process: first, we estimate parameters, second (below) we filter the data.
rule variant_recalibrator:
    input:
        ref=config["data"]["reference-genome"],
        vcf="filtered/all.select-{vartype}.vcf.gz",

        # Resources have to be given as named input files, we unpack the dict to get them.
        # We also request the tbi index files, so that users get a nice error message if missing.
        **get_variant_recalibrator_resource_files(),
        tbi=get_variant_recalibrator_resource_file_tbis(),

    output:
        # Ouput file needs to be called vcf for the wrapper, but is in fact a fake vcf
        # that actually contains the recal information.
        vcf="filtered/all.{vartype}.vqsr-recal.vcf.gz",
        tranches="filtered/all.{vartype}.vqsr-recal.tranches"
        # We also might produce a plot PDF about the trances - but only for the SNPs,
        # not for the INDELs, so we do not specify it here for simplicity...
        # The avid user will find it in the `filtered` directory either way.
    params:
        # set mode, must be either SNP, INDEL or BOTH
        mode="{vartype}",

        # Resource parameter definition. Key must match named input files from above.
        # Snakemake/Python makes this a bit difficult with the config file...
        resources = get_variant_recalibrator_resources(),

        # which fields to use with -an (see VariantRecalibrator docs)
        annotation=config["params"]["variantfiltration-vqsr"]["annotation"],

        # Extras
        extra=get_variant_recalibrator_extra,
        java_opts=config["params"]["variantfiltration-vqsr"]["java-variantrecalibrator"],
    log:
        "logs/gatk/variantrecalibrator/{vartype}.log"
    benchmark:
        "benchmarks/gatk/variantrecalibrator/{vartype}.bench.log"
    # Group deactivated, so that it can run in parallel for SNP and INDEL
    # group:
    #     "filtering"
    resources:
        # optional specification of memory usage of the JVM that snakemake will respect with global
        # resource restrictions
        # (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
        # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
        # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
        # mem_mb=1024
    conda:
        # We overwrite the original yaml, as it does not contain the R specifications,
        # which we want in order to also plot the trances file for SNPs... because why not.
        "../envs/gatk.yaml"
    wrapper:
        "0.85.0/bio/gatk/variantrecalibrator"

rule apply_vqsr:
    input:
        ref=config["data"]["reference-genome"],
        vcf="filtered/all.select-{vartype}.vcf.gz",
        recal="filtered/all.{vartype}.vqsr-recal.vcf.gz",
        tranches="filtered/all.{vartype}.vqsr-recal.tranches"
    output:
        vcf=(
            "filtered/all.{vartype}.recalibrated.vcf.gz"
            if config["settings"]["keep-intermediate"]["filtering"]
            else temp("filtered/all.{vartype}.recalibrated.vcf.gz")
        )
    log:
        "logs/gatk/applyvqsr/{vartype}.log"
    benchmark:
        "benchmarks/gatk/applyvqsr/{vartype}.bench.log"
    params:
        # set mode, must be either SNP, INDEL or BOTH
        mode="{vartype}",
        extra=get_apply_vqsr_extra,
        java_opts=config["params"]["variantfiltration-vqsr"]["java-applyvqsr"]
    resources:
        # mem_mb=50
    conda:
        # We overwrite the original yaml, as this wrapper here (version 0.85.0) and the one above
        # for the variantrecalibrator (also 0.85.0) use different GATK versions originally...
        # Ah if only things were consistent... so let's make them use the same version here.
        "../envs/gatk.yaml"
    wrapper:
        "0.85.0/bio/gatk/applyvqsr"

# =================================================================================================
#     Merge Filtered Variants
# =================================================================================================

rule merge_calls:
    input:
        vcf=expand(
            "filtered/all.{vartype}.{filtertype}.vcf.gz",
            vartype=["SNP", "INDEL"],
            filtertype="recalibrated"
            if config["settings"].get("vqsr", False)
            else "hardfiltered"
        )
    output:
        vcf="filtered/all.vcf.gz"
        # vcf=protected("filtered/all.vcf.gz")
    log:
        "logs/picard/merge-filtered.log"
    benchmark:
        "benchmarks/picard/merge-filtered.bench.log"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"
