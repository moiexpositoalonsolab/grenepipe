import json
import platform

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
            "filtered/all.{vartype}.selected.vcf.gz"
            if config["settings"]["keep-intermediate"]["filtering"]
            else temp("filtered/all.{vartype}.selected.vcf.gz")
        ),
        done=touch("filtered/all.{vartype}.selected.done"),
    params:
        extra="--select-type-to-include {vartype}",
    log:
        "logs/gatk/selectvariants/{vartype}.log",
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

# Switch to the chosen filtering tool
filter_variants_good = False
if config["settings"]["filter-variants"] == "gatk-variantfiltration":

    # Use `gatk-variantfiltration`
    filter_variants_good = True

    include: "filtering-gatk-variantfiltration.smk"

elif config["settings"]["filter-variants"] == "gatk-vqsr":

    # Use `gatk-vqsr`
    filter_variants_good = True

    include: "filtering-gatk-vqsr.smk"

elif config["settings"]["filter-variants"] == "bcftools-filter":

    # Use `bcftools-filter`
    filter_variants_good = True

    include: "filtering-bcftools-filter.smk"


if config["settings"]["filter-variants"] == "none":

    # Nothing to include
    filter_variants_good = True

# Slightly awkward setup of the conditions here, due to a bug in snakefmt,
# see https://github.com/snakemake/snakefmt/issues/239#issuecomment-2223276169
if not filter_variants_good:
    raise Exception("Unknown filter-variants: " + config["settings"]["filter-variants"])

# =================================================================================================
#     Merge Filtered Variants
# =================================================================================================


rule merge_calls:
    input:
        # We use different naming for the VQSR intermediate files,
        # to make it a bit more understandable to the user which file is which...
        # not sure if that helps, or is more confusing in the end :-O
        vcf=expand(
            "filtered/all.{vartype}.{filtertype}.vcf.gz",
            vartype=["SNP", "INDEL"],
            filtertype=(
                "recalibrated"
                if config["settings"]["filter-variants"] == "gatk-vqsr"
                else "filtered"
            ),
        ),
    output:
        vcf="filtered/all.vcf.gz",
        # vcf=protected("filtered/all.vcf.gz")
        done=touch("filtered/all.done"),
    params:
        # See duplicates-picard.smk for the reason whe need this on MacOS.
        extra=(
            " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true" if platform.system() == "Darwin" else ""
        ),
    log:
        "logs/picard/merge-filtered.log",
    benchmark:
        "benchmarks/picard/merge-filtered.bench.log"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"
