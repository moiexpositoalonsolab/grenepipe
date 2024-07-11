# =================================================================================================
#     bcftools filter
# =================================================================================================

# Need a function, in order to use wildcards within params dicts...
def get_filter(wildcards):
    return config["params"]["bcftools-filter"][wildcards.vartype]

rule bcftools_filter_calls:
    input:
        vcf="filtered/all.{vartype}.selected.vcf.gz"
    output:
        vcf=(
            "filtered/all.{vartype}.filtered.vcf.gz"
            if config["settings"]["keep-intermediate"]["filtering"]
            else temp("filtered/all.{vartype}.filtered.vcf.gz")
        ),
        done=touch("filtered/all.{vartype}.filtered.done")
    params:
        tool=config["params"]["bcftools-filter"]["tool"],
        filters=get_filter,
        extra=config["params"]["bcftools-filter"]["extra"]
    log:
        "logs/bcftools/filter-{vartype}.log"
    benchmark:
        "benchmarks/bcftools/filter-{vartype}.bench.log"
    group:
        "filtering"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "("
        "bcftools {params.tool} {params.filters} {params.extra} -O z -o {output.vcf} {input.vcf} ; "
        "bcftools index --tbi {output.vcf}"
        ") &> {log}"
