# =================================================================================================
#     Hard Filtering
# =================================================================================================

# Need a function, in order to use wildcards within params dicts...
def get_filter(wildcards):
    return {
        wildcards.vartype + "-hard-filter":
        config["params"]["gatk-variantfiltration"][wildcards.vartype]
    }

# Simple hard filter, used by default
rule gatk_hard_filter_calls:
    input:
        ref=config["data"]["reference-genome"],
        vcf="filtered/all.{vartype}.selected.vcf.gz",
        refdict=genome_dict()
    output:
        vcf=(
            "filtered/all.{vartype}.filtered.vcf.gz"
            if config["settings"]["keep-intermediate"]["filtering"]
            else temp("filtered/all.{vartype}.filtered.vcf.gz")
        ),
        done=touch("filtered/all.{vartype}.filtered.done")
    params:
        filters=get_filter,
        extra=config["params"]["gatk-variantfiltration"]["extra"],
        java_opts=config["params"]["gatk-variantfiltration"]["java-opts"]
    log:
        "logs/gatk/variantfiltration/{vartype}.log"
    benchmark:
        "benchmarks/gatk/variantfiltration/{vartype}.bench.log"
    group:
        "filtering"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "0.85.0/bio/gatk/variantfiltration"
