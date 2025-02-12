# =================================================================================================
#     Base Quality Score Recalibration
# =================================================================================================

# If recal base quals is set, we need a few extra checks to make sure that this works.
# This gives nicer user output in cases of these issues than relying on the tool-internal errors.
# The platform information is provided for GATK BaseRecalibrator via the bam RG tags,
# meaning that the mapping tool needs to know them and add them to the bam files.
if config["settings"]["recalibrate-base-qualities"]:
    if not config["data"]["known-variants"]:
        raise Exception(
            "Setting recalibrate-base-qualities can only be activated if a known-variants "
            "file is also provided, as GATK BaseRecalibrator needs this."
        )
    if "platform" not in config["global"]["samples"] and not config["params"]["gatk"]["platform"]:
        raise Exception(
            "Setting recalibrate-base-qualities can only be activated if either a `platform` column "
            "is provided in the input samples table ("
            + config["global"]["samples"]
            + ") or if a platform for all samples is given in the `params: gatk: platform` setting "
            "in the config file."
        )
    if "platform" in config["global"]["samples"] and config["params"]["gatk"]["platform"]:
        logger.warning(
            "The samples table contains a column `platform` for each sample, "
            "but the `params: gatk: platform` setting is also provided. "
            "We will use the table, as this is more specific."
        )
    if config["params"]["gatk"]["platform"]:
        platforms = set(config["params"]["gatk"]["platform"])
    if "platform" in config["global"]["samples"]:
        platforms = set(config["global"]["samples"]["platform"])
    for p in ["CAPILLARY", "LS454", "ILLUMINA", "SOLID", "HELICOS", "IONTORRENT", "ONT", "PACBIO"]:
        if p in platforms:
            platforms.remove(p)
    if platforms:
        logger.warning(
            "Provided sequencing platforms (in the samples table `platform` column, or in the "
            "`params: gatk: platform` setting) contains values that might not be recognized by "
            "GATK BaseRecalibrator, and hence might lead to errors: " + str(platforms)
        )


def get_recal_input(ext=""):
    # case 1: no duplicate removal
    f = "mapping/merged/{sample}.bam"

    # case 2: filtering via samtools view
    if config["settings"]["filter-mapped-reads"]:
        f = "mapping/filtered/{sample}.bam"

    # case 3: clipping reads with BamUtil
    if config["settings"]["clip-read-overlaps"]:
        f = "mapping/clipped/{sample}.bam"

    # case 4: remove duplicates
    if config["settings"]["remove-duplicates"]:
        f = "mapping/dedup/{sample}.bam"

    if ext == ".bai":
        if config["settings"].get("restrict-regions"):
            # case 5: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 6: no index needed
            return []
    else:
        return f + ext


# Helper function needed below. Could probably be replaced by intervals.
def get_gatk_regions_param(regions=config["settings"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        # Not used at the moment, as we deleted this config setting.
        # padding = config["settings"].get("region-padding")
        # if padding:
        #     params += "--interval-padding {}".format(padding)
        return params
    return default


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(ext=".bai"),
        done=get_recal_input(ext=".done"),
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        refdict=genome_dict(),
        known=config["data"]["known-variants"],
    output:
        bam=(
            "mapping/recal/{sample}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapping/recal/{sample}.bam")
        ),
        # bam=protected("mapping/recal/{sample}.bam")
        done=touch("mapping/recal/{sample}.bam.done"),
    params:
        extra=get_gatk_regions_param() + " " + config["params"]["gatk"]["BaseRecalibrator"],
    log:
        "logs/mapping/gatk-bqsr/{sample}.log",
    benchmark:
        "benchmarks/mapping/gatk-bqsr/{sample}.log"
    group:
        "mapping_extra"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "0.51.3/bio/gatk/baserecalibrator"
