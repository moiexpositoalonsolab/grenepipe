import json

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
    rf = ordered_dict_to_dict(config["params"]["gatk-vqsr"]["resource-files"])
    if any(x in rf.keys() for x in ["ref", "vcf", "tbi"]):
        raise Exception("Cannot use resource named 'ref', 'vcf', or 'tbi' with our GATK VQSR.")
    return rf


# We also add the tbi files as dependencies, so that their existence is checked.
def get_variant_recalibrator_resource_file_tbis():
    rf = get_variant_recalibrator_resource_files()
    tbi = []
    for k in rf.keys():
        tbi.append(rf[k] + ".tbi")
    return tbi


def get_variant_recalibrator_resources():
    res = ordered_dict_to_dict(config["params"]["gatk-vqsr"]["resources"])
    for rn, rv in res.items():
        if isinstance(rv, list):
            nv = {}
            for erv in rv:
                for k, v in erv.items():
                    nv[k] = v
            res[rn] = nv
    return res


def get_variant_recalibrator_extra(wildcards):
    # Need to do this in a function, as wildcard replacement otherwise does not work within
    # config file dict retrieval. We additionally set the Rscript file, so that we get some plots.
    vartype = wildcards.vartype
    return (
        config["params"]["gatk-vqsr"]["variantrecalibrator-extra-" + vartype]
        + " --rscript-file calling/filtered/all."
        + vartype
        + ".vqsr-recal.plots.R"
    )


def get_apply_vqsr_extra(wildcards):
    # Same as above, need a function here for wildcard replacement
    return config["params"]["gatk-vqsr"]["applyvqsr-extra-" + wildcards.vartype]


# Use the GATK VQSR machine learning based recalibration of quality scores instead of hard filtering.
# This is a two step process: first, we estimate parameters, second (below) we filter the data.
rule gatk_variant_recalibrator:
    input:
        **get_variant_recalibrator_resource_files(),
        ref=config["data"]["reference-genome"],
        vcf="calling/filtered/all.{vartype}.selected.vcf.gz",
        done="calling/filtered/all.{vartype}.selected.vcf.gz.done",
        refdict=genome_dict(),
        # Resources have to be given as named input files, we unpack the dict to get them.
        # We also request the tbi index files, so that users get a nice error message if missing.
        tbi=get_variant_recalibrator_resource_file_tbis(),
    output:
        # Ouput file needs to be called vcf for the wrapper, but is in fact a fake vcf
        # that actually contains the recal information.
        vcf="calling/filtered/all.{vartype}.vqsr-recal.vcf.gz",
        tranches="calling/filtered/all.{vartype}.vqsr-recal.tranches",
        # We also might produce a plot PDF about the trances - but only for the SNPs,
        # not for the INDELs, so we do not specify it here for simplicity...
        # The avid user will find it in the `filtered` directory either way.
        done=touch("calling/filtered/all.{vartype}.vqsr-recal.vcf.gz.done"),
    params:
        # set mode, must be either SNP, INDEL or BOTH
        mode="{vartype}",
        # Resource parameter definition. Key must match named input files from above.
        # Snakemake/Python makes this a bit difficult with the config file...
        resources=get_variant_recalibrator_resources(),
        # which fields to use with -an (see VariantRecalibrator docs)
        annotation=config["params"]["gatk-vqsr"]["annotation"],
        # Extras
        extra=get_variant_recalibrator_extra,
        java_opts=config["params"]["gatk-vqsr"]["variantrecalibrator-java-opts"],
    resources:
        mem_mb=config["params"]["gatk-vqsr"].get("variantrecalibrator-mem-mb", 1024),
    log:
        "logs/calling/gatk-variantrecalibrator/{vartype}.log",
    benchmark:
        "benchmarks/calling/gatk-variantrecalibrator/{vartype}.log"
    # Group deactivated, so that it can run in parallel for SNP and INDEL
    # group:
    #     "filtering"
    # resources:
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


rule gatk_apply_vqsr:
    input:
        ref=config["data"]["reference-genome"],
        refdict=genome_dict(),
        vcf="calling/filtered/all.{vartype}.selected.vcf.gz",
        vcf_done="calling/filtered/all.{vartype}.selected.vcf.gz.done",
        recal="calling/filtered/all.{vartype}.vqsr-recal.vcf.gz",
        tranches="calling/filtered/all.{vartype}.vqsr-recal.tranches",
        recal_done="calling/filtered/all.{vartype}.vqsr-recal.vcf.gz.done",
    output:
        vcf=(
            "calling/filtered/all.{vartype}.recalibrated.vcf.gz"
            if config["settings"]["keep-intermediate"]["filtering"]
            else temp("calling/filtered/all.{vartype}.recalibrated.vcf.gz")
        ),
        done=touch("calling/filtered/all.{vartype}.recalibrated.vcf.gz.done"),
    log:
        "logs/calling/gatk-applyvqsr/{vartype}.log",
    benchmark:
        "benchmarks/calling/gatk-applyvqsr/{vartype}.log"
    params:
        # set mode, must be either SNP, INDEL or BOTH
        mode="{vartype}",
        extra=get_apply_vqsr_extra,
        java_opts=config["params"]["gatk-vqsr"]["applyvqsr-java-opts"],
    resources:
        mem_mb=config["params"]["gatk-vqsr"].get("applyvqsr-mem-mb", 1024),
    conda:
        # We overwrite the original yaml, as this wrapper here (version 0.85.0) and the one above
        # for the variantrecalibrator (also 0.85.0) use different GATK versions originally...
        # Ah if only things were consistent... so let's make them use the same version here.
        "../envs/gatk.yaml"
    wrapper:
        "0.85.0/bio/gatk/applyvqsr"
