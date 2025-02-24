# =================================================================================================
#     Hard Filtering
# =================================================================================================


# Need a function, in order to use wildcards within params dicts...
def get_filter(wildcards):
    return {
        wildcards.vartype
        + "-hard-filter": config["params"]["gatk-variantfiltration"][wildcards.vartype]
    }


# Simple hard filter, used by default
rule gatk_hard_filter_calls:
    input:
        ref=config["data"]["reference-genome"],
        vcf="calling/filtered/all.{vartype}.selected.vcf.gz",
        don="calling/filtered/all.{vartype}.selected.vcf.gz.done",
        refdict=genome_dict(),
    output:
        vcf=(
            "calling/filtered/all.{vartype}.filtered.vcf.gz"
            if config["settings"]["keep-intermediate"]["filtering"]
            else temp("calling/filtered/all.{vartype}.filtered.vcf.gz")
        ),
        done=touch("calling/filtered/all.{vartype}.filtered.vcf.gz.done"),
    params:
        filters=get_filter,
        extra=config["params"]["gatk-variantfiltration"]["extra"],
        java_opts=config["params"]["gatk-variantfiltration"]["java-opts"],
    resources:
        mem_mb=config["params"]["gatk-variantfiltration"].get("mem-mb", 1024),
    log:
        "logs/calling/gatk-variantfiltration/{vartype}.log",
    benchmark:
        "benchmarks/calling/gatk-variantfiltration/{vartype}.log"
    group:
        "filtering"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/variantfiltration"
