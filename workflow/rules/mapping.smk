# =================================================================================================
#     Read Group and Helper Functions
# =================================================================================================


# Generic rule for all bai files.
# Bit weird as it produces log files with nested paths, but that's okay for now.
# To avoid having absolute paths as part of the log file output though,
# which can happen for the absolute paths provided by the mappings table,
# we add an underscore, just to silence the snakemake warning...
# Super ugly, but log file paths do not work with lambdas or functions in snakemake,
# so there is no easy solution. And introducing separate bam index functions
# for every case where we need those is stupid as well.
# So we live with the ugly understore and neste paths... It's just a log file,
# after all. Still, ugly due to snakemake being weird.
rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        # "logs/mapping/samtools-index/{prefix}.log",
        "logs/mapping/samtools-index/_" + "{prefix}.log",
    group:
        "mapping_extra"
    wrapper:
        "0.51.3/bio/samtools/index"


# We here get the read group tags list that is used per sample for the reads in the bam file.
# This differs a bit depending on whether the `platform` field is properly set in the samples table.
# We do not construct a full string here, as mapping tools take this information in different ways.
def get_read_group_tags(wildcards):
    # We need the @RG read group tags, including `ID` and `SM`, as downstream tools use these.
    # Potentially, we also set the PL platform. Also add the config extra settings.
    res = ["ID:" + wildcards.sample + "-" + wildcards.unit, "SM:" + wildcards.sample]
    # TODO Add LD? LB? field as well for the unit?! http://www.htslib.org/workflow/

    # Add platform information, if available, giving precedence to the table over the config.
    pl = ""
    if "platform" in config["params"]["gatk"] and config["params"]["gatk"]["platform"]:
        pl = config["params"]["gatk"]["platform"]
    if "platform" in config["global"]["samples"]:
        s = (
            config["global"]["samples"]
            .loc[(wildcards.sample, wildcards.unit), ["platform"]]
            .dropna()
        )
        # We catch the case that the platform column is present, but empty for a given row.
        # This would lead to an error, see https://stackoverflow.com/q/610883/4184258
        try:
            pl = s.platform
        except:
            pl = ""
    if pl:
        res.append("PL:" + pl)
    return res

    # For reference: https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/blob/cffa77a9634801152971448bd41d0687cf765723/workflow/rules/common.smk#L60
    # return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
    #     sample=wildcards.sample,
    #     platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
    # )


# Helper function needed by some of the GATK tools, here and in the calling.
def get_gatk_regions_param(regions=config["settings"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        # Not used at the moment, as we deleted this config setting.
        # padding = config["settings"].get("region-padding")
        # if padding:
        #     params += "--interval-padding {}".format(padding)
        return params
    return default


# =================================================================================================
#     Read Mapping
# =================================================================================================

# Switch to the chosen mapper
mapping_tool_good = False
if config["settings"]["mapping-tool"] == "bwaaln":

    # Use `bwa aln`
    include: "mapping-bwa-aln.smk"

    mapping_tool_good = True

elif config["settings"]["mapping-tool"] == "bwamem":

    # Use `bwa mem`
    include: "mapping-bwa-mem.smk"

    mapping_tool_good = True

elif config["settings"]["mapping-tool"] == "bwamem2":

    # Use `bwa mem2`
    include: "mapping-bwa-mem2.smk"

    mapping_tool_good = True

elif config["settings"]["mapping-tool"] == "bowtie2":

    # Use `bowtie2`
    include: "mapping-bowtie2.smk"

    mapping_tool_good = True


# Weird syntax instead of just using an else branch here, due to
# https://github.com/snakemake/snakefmt/issues/239#issuecomment-2223276169
if not mapping_tool_good:
    raise Exception("Unknown mapping-tool: " + config["settings"]["mapping-tool"])

# =================================================================================================
#     Merge per-sample bam files
# =================================================================================================


# We need a helper function to expand based on wildcards.
# The file names used here is what all the above mappers are expected to produce.
def get_sorted_sample_bams(wildcards):
    return expand("mapping/sorted/{{sample}}-{unit}.bam", unit=get_sample_units(wildcards.sample))


# This is where all units are merged together.
# We changed this behaviour. grenepipe v0.11.1 and before did _all_ the steps in this file
# individually per unit. After that, we changed to merge here, and then run every subsequent
# step on the samples with all their units merged. This makes the processing faster, and ensures
# that all tools properly understand that these files are supposed to be a single sample each,
# instead of potentially accidentally considering the different units as individual samples.
rule merge_sample_unit_bams:
    input:
        get_sorted_sample_bams,
    output:
        "mapping/merged/{sample}.bam",
        touch("mapping/merged/{sample}.done"),
    params:
        extra=config["params"]["samtools"]["merge"],
    threads: config["params"]["samtools"]["merge-threads"]
    log:
        "logs/mapping/samtools-merge/{sample}.log",
    wrapper:
        "v3.13.6/bio/samtools/merge"


# =================================================================================================
#     Filtering and Clipping Mapped Reads
# =================================================================================================


rule filter_mapped_reads:
    input:
        "mapping/merged/{sample}.bam",
    output:
        (
            "mapping/filtered/{sample}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapping/filtered/{sample}.bam")
        ),
        touch("mapping/filtered/{sample}.done"),
    params:
        extra=config["params"]["samtools"]["view"] + " -b",
    log:
        "logs/mapping/samtools-view/{sample}.log",
    conda:
        # Need our own env again, because of conflicting numpy and pandas version...
        "../envs/samtools.yaml"
    wrapper:
        "0.85.0/bio/samtools/view"


rule clip_read_overlaps:
    input:
        # Either get the mapped reads directly, or if we do filtering before, take that instead.
        "mapping/filtered/{sample}.bam"
        if (config["settings"]["filter-mapped-reads"])
        else ("mapping/merged/{sample}.bam"),
    output:
        (
            "mapping/clipped/{sample}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapping/clipped/{sample}.bam")
        ),
        touch("mapping/clipped/{sample}.done"),
    params:
        extra=config["params"]["bamutil"]["extra"],
    log:
        "logs/mapping/bamutil/{sample}.log",
    benchmark:
        "benchmarks/mapping/bamutil/{sample}.log"
    conda:
        "../envs/bamutil.yaml"
    shell:
        "bam clipOverlap --in {input[0]} --out {output[0]} {params.extra} &> {log}"


def get_mapped_reads(wildcards):
    """Get mapped reads of given sample (with merged units),
    either unfiltered, filtered, and/or clipped."""

    # Default case: just the mapped and sorted reads.
    result = "mapping/merged/{sample}.bam".format(**wildcards)

    # If we want filtering, do that instead.
    if config["settings"]["filter-mapped-reads"]:
        result = "mapping/filtered/{sample}.bam".format(**wildcards)

    # If we want clipping, do that. The clipping rule above also might already use the
    # filtered reads, so that is taking into account as well.
    if config["settings"]["clip-read-overlaps"]:
        result = "mapping/clipped/{sample}.bam".format(**wildcards)

    return result


# We need to get a bit dirty here, as dedup names output files based on input file names,
# so that depending on what kind of input (mapped, filtered, clipped) we give it,
# we get differently named output files, which does not work well with snakemake.
# So we need to hardcode this here, it seems...
# We keep this function here, so that if we ever add an additinal step here, we will (hopefully)
# notice, and adapt this function as well.
# Udate: We now restructured the directories, so that the file names are consistent,
# and instead have differently named parent directories for the individual steps, see above.
# So now it seems we do not need this funciton here any more. Still kept for reference.
# def get_mapped_read_infix():
#     # Same logic as the function above.
#     result = "merged"
#     if config["settings"]["filter-mapped-reads"]:
#         result = "filtered"
#     if config["settings"]["clip-read-overlaps"]:
#         result = "clipped"
#     return result


# =================================================================================================
#     Mark Duplicates
# =================================================================================================

# Switch to the chosen duplicate marker tool
duplicates_tool_good = False
if config["settings"]["duplicates-tool"] == "picard":

    # Bioinformatics tools are messed up...
    if config["settings"]["mapping-tool"] == "bowtie2":
        raise Exception(
            "Cannot combine mapping-tool bowtie2 with duplicates-tool picard, "
            "because those two have not yet learned to work with each others file formats. "
            "In particular, bowtie2 produces a @PG header line in its output bam file "
            "that picard cannot parse. Unclear who is to blame. "
            "If you need this combination of tools, please submit an issue to "
            "https://github.com/lczech/grenepipe/issues and we will see what we can do."
        )
        # See also https://github.com/broadinstitute/picard/issues/1139
        # and https://github.com/samtools/htsjdk/issues/677
        # for previous encounters of this issue. Seems not solved yet, so we disallow this for now.

    # Use `picard`
    include: "duplicates-picard.smk"

    duplicates_tool_good = True

elif config["settings"]["duplicates-tool"] == "dedup":

    # Use `dedup`
    include: "duplicates-dedup.smk"

    duplicates_tool_good = True


# Again, fix for https://github.com/snakemake/snakefmt/issues/239#issuecomment-2223276169
if not duplicates_tool_good:
    raise Exception("Unknown duplicates-tool: " + config["settings"]["duplicates-tool"])

# =================================================================================================
#     Base Quality Score Recalibration
# =================================================================================================

if config["settings"]["recalibrate-base-qualities"]:

    include: "mapping-recalibrate.smk"


# =================================================================================================
#     Final Mapping Result
# =================================================================================================


# Helper function for the case that a `mappings-table` is provided, in which case we do not run
# any mapping ourselves, and skipp all of the above. Instead, we then simply want to use the
# sample bam file from the table here.
def get_bam_from_mappings_table(sample):
    return config["global"]["samples"].loc[sample, ["bam"]].dropna()


# At this point, we have several choices of which files we want to hand down to the next
# pipleine step. We offer a function so that downstream does not have to deal with this.
# The way this works is as follows: Each optional step that could be applied
# overwrites the previous setting, until we have the final file that we want to use downstream.
# As each of these optional steps itself also takes care of using previous optional steps
# (see above), what we get here is the file name of the last step that we want to apply.


def get_mapping_result(sample, bai=False):
    # Special case: we are using the mappings table of bam files,
    # instead of any of the ones produced with the rules here.
    if "mappings-table" in config["data"] and config["data"]["mappings-table"]:
        f = get_bam_from_mappings_table(sample)

    else:
        # All other cases: pick the "final" bam file following the settings
        # that is the one to be used downstream from here.

        # case 1: no duplicate removal
        f = "mapping/merged/{sample}.bam".format(sample=sample)

        # case 2: filtering via samtools view
        if config["settings"]["filter-mapped-reads"]:
            f = "mapping/filtered/{sample}.bam".format(sample=sample)

        # case 3: clipping reads with BamUtil
        if config["settings"]["clip-read-overlaps"]:
            f = "mapping/clipped/{sample}.bam".format(sample=sample)

        # case 4: remove duplicates
        if config["settings"]["remove-duplicates"]:
            f = "mapping/dedup/{sample}.bam".format(sample=sample)

        # case 5: recalibrate base qualities
        if config["settings"]["recalibrate-base-qualities"]:
            f = "mapping/recal/{sample}.bam".format(sample=sample)

    # Additionally, this function is run for getting bai files as well
    if bai:
        f += ".bai"

    return f


# Return the bam file for a given sample.
# Get all aligned reads of given sample, with all its units merged.
# The function automatically gets which of the mapping results to use, depending on the config
# setting (whether to remove duplicates, and whether to recalibrate the base qualities),
# by using the get_mapping_result function, that gives the respective files depending on the config.
def get_sample_bams(sample):
    return get_mapping_result(sample)


# Return the bai file(s) for a given sample
def get_sample_bais(sample):
    return get_mapping_result(sample, True)


# Return the bam file(s) for a sample, given a wildcard param from a rule.
def get_sample_bams_wildcards(wildcards):
    return get_sample_bams(wildcards.sample)


# Return the bai file(s) for a sample, given a wildcard param from a rule.
def get_sample_bais_wildcards(wildcards):
    return get_sample_bais(wildcards.sample)


# Return the bam file(s) for all samples
def get_all_bams():
    # Make a list of all bams in the order as the samples list.
    res = list()
    for sample in config["global"]["sample-names"]:
        res.append(get_mapping_result(sample))
    return res


# Return the bai file(s) for all samples
def get_all_bais():
    res = list()
    for sample in config["global"]["sample-names"]:
        res.append(get_mapping_result(sample, True))
    return res


# =================================================================================================
#     All bams, but not SNP calling
# =================================================================================================


# This alternative target rule executes all steps up to th mapping, and yields the final bam
# files that would otherwise be used for variant calling in the downstream process.
# That is, depending on the config, these are the sorted+merged, filtered, remove duplicates, or
# recalibrated base qualities bam files.

# Additionally, we create a "final" mapping result dir, for user convenience,
# which simply contains symlinks to whatever the final result is.
# That way, if a user is looking for the final result of the mapping step,
# they just have to consult this one.


rule all_bams:
    input:
        bams=get_all_bams(),
        qc="qc/multiqc.html",
    output:
        bams=expand("mapping/final/{sample}.bam", sample=config["global"]["sample-names"]),
        done=touch("mapping/final.done"),
    run:
        for bam in input.bams:
            sample = os.path.basename(bam)
            output_bam = f"mapping/final/{sample}"
            bam_abs_path = os.path.abspath(bam)
            if not os.path.exists(output_bam):
                shell(f"ln -s {bam_abs_path} {output_bam}")


# The `all_bams` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules:
    all_bams,
