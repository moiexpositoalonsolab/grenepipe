import platform

# =================================================================================================
#     Common Functions used for bcftools calling
# =================================================================================================

if config["data"]["known-variants"]:
    raise Exception(
        "Calling tool 'bcftools' cannot be used with the option 'known-variants', "
        "as it does not support to call based on known variants."
    )


def get_mpileup_params(wildcards, input):
    # Start with the user-specified params from the config file
    params = config["params"]["bcftools"]["mpileup"]

    # input.regions is either the restrict-regions or small contigs bed file, or empty.
    # If given, we need to set a param of mpileup to use the file.
    # If not given, we set a param that uses the contig string instead.
    if input.regions:
        params += " --regions-file {}".format(input.regions)
    else:
        params += " --regions {}".format(wildcards.contig)

    # Lastly, we need some extra annotations, so that our downstream stats rules can do their job.
    # That was a hard debugging session to figure this one out... If only bioinformatics tools
    # had better error messages!
    params += " --annotate FORMAT/AD,FORMAT/DP"
    return params


# =================================================================================================
#     Variant Calling
# =================================================================================================

# Switch to the chosen caller mode
bcftools_mode_good = False
if config["params"]["bcftools"]["mode"] == "combined":

    bcftools_mode_good = True

    include: "calling-bcftools-combined.smk"

elif config["params"]["bcftools"]["mode"] == "individual":

    bcftools_mode_good = True

    include: "calling-bcftools-individual.smk"


# Stupid workaround for https://github.com/snakemake/snakefmt/issues/239
if not bcftools_mode_good:
    raise Exception("Unknown bcftools mode: " + config["params"]["bcftools"]["mode"])

# =================================================================================================
#     Sorting Calls
# =================================================================================================

# When using small contigs, we further need to sort the output, as
# this won't be done for us. Let's only do this extra work though if needed.
if config["settings"].get("contig-group-size"):

    rule sort_variants:
        input:
            vcf="calling/merged/merged-all.vcf.gz",
            done="calling/merged/merged-all.vcf.gz.done",
            refdict=genome_dict(),
        output:
            vcf="calling/genotyped-all.vcf.gz",
            done=touch("calling/genotyped-all.vcf.gz.done"),
        params:
            # See duplicates-picard.smk for the reason whe need this on MacOS.
            extra=(
                " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
                if platform.system() == "Darwin"
                else ""
            ),
            java_opts=config["params"]["picard"]["SortVcf-java-opts"],
        log:
            "logs/calling/picard/sort-genotyped.log",
        benchmark:
            "benchmarks/calling/genotyped/picard/sort-genotyped.log"
        conda:
            "../envs/picard.yaml"
        shell:
            # Weird new picard syntax...
            "picard SortVcf "
            "{params.java_opts} "
            "INPUT={input.vcf} "
            "OUTPUT={output.vcf} "
            "SEQUENCE_DICTIONARY={input.refdict} "
            "{params.extra} "
            "> {log} 2>&1"
