from itertools import chain

# =================================================================================================
#     Read Group
# =================================================================================================

# We here get the read group tags list that is used per sample for the reads in the bam file.
# This differs a bit depending on whether the `platform` field is properly set in the samples table.
# We do not construct a full string here, as mapping tools take this information in different ways.
def get_read_group_tags( wildcards ):
    # We need the @RG read group tags, including `ID` and `SM`, as downstream tools use these.
    # Potentially, we also set the PL platform. Also add the config extra settings.
    res = [ "ID:" + wildcards.sample, "SM:" + wildcards.sample ]
    # TODO Add LD field as well for the unit?! http://www.htslib.org/workflow/

    # Add platform information, if available, giving precedence to the table over the config.
    pl = ""
    if 'platform' in config["params"]["gatk"] and config["params"]["gatk"]["platform"]:
        pl = config["params"]["gatk"]["platform"]
    if 'platform' in config["global"]["samples"]:
        s = config["global"]["samples"].loc[(wildcards.sample, wildcards.unit), ["platform"]].dropna()
        pl = s.platform
    if pl:
        res.append( "PL:" + pl )
    return res

    # For reference: https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/blob/cffa77a9634801152971448bd41d0687cf765723/workflow/rules/common.smk#L60
    # return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
    #     sample=wildcards.sample,
    #     platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
    # )

# =================================================================================================
#     Read Mapping
# =================================================================================================

# Switch to the chosen mapper
if config["settings"]["mapping-tool"] == "bwaaln":

    # Use `bwa aln`
    include: "mapping-bwa-aln.smk"

elif config["settings"]["mapping-tool"] == "bwamem":

    # Use `bwa mem`
    include: "mapping-bwa-mem.smk"

elif config["settings"]["mapping-tool"] == "bwamem2":

    # Use `bwa mem2`
    include: "mapping-bwa-mem2.smk"

elif config["settings"]["mapping-tool"] == "bowtie2":

    # Use `bowtie2`
    include: "mapping-bowtie2.smk"

else:
    raise Exception("Unknown mapping-tool: " + config["settings"]["mapping-tool"])

# =================================================================================================
#     Filtering and Clipping Mapped Reads
# =================================================================================================

rule filter_mapped_reads:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        (
            "mapped/{sample}-{unit}.filtered.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapped/{sample}-{unit}.filtered.bam")
        ),
        touch("mapped/{sample}-{unit}.filtered.done")
    params:
        extra=config["params"]["samtools"]["view"] + " -b"
    conda:
        # Need our own env again, because of conflicting numpy and pandas version...
        "../envs/samtools.yaml"
    wrapper:
        "0.85.0/bio/samtools/view"

rule clip_read_overlaps:
    input:
        # Either get the mapped reads directly, or if we do filtering before, take that instead.
        "mapped/{sample}-{unit}.filtered.bam" if (
            config["settings"]["filter-mapped-reads"]
        ) else (
            "mapped/{sample}-{unit}.sorted.bam"
        )
    output:
        (
            "mapped/{sample}-{unit}.clipped.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapped/{sample}-{unit}.clipped.bam")
        ),
        touch("mapped/{sample}-{unit}.clipped.done")
    params:
        extra=config["params"]["bamutil"]["extra"]
    log:
        "logs/bamutil/{sample}-{unit}.log"
    benchmark:
        "benchmarks/bamutil/{sample}-{unit}.bench.log"
    conda:
        "../envs/bamutil.yaml"
    shell:
        "bam clipOverlap --in {input[0]} --out {output[0]} {params.extra} &> {log}"

def get_mapped_reads(wildcards):
    """Get mapped reads of given sample-unit, either unfiltered, filtered, and/or clipped."""

    # Default case: just the mapped and sorted reads.
    result = "mapped/{sample}-{unit}.sorted.bam".format(**wildcards)

    # If we want filtering, do that instead.
    if config["settings"]["filter-mapped-reads"]:
        result = "mapped/{sample}-{unit}.filtered.bam".format(**wildcards)

    # If we want clipping, do that. The clipping rule above also might already use the
    # filtered reads, so that is taking into account as well.
    if config["settings"]["clip-read-overlaps"]:
        result = "mapped/{sample}-{unit}.clipped.bam".format(**wildcards)

    return result

# We need to get a bit dirty here, as dedup names output files based on input file names,
# so that depending on what kind of input (mapped, filtered, clipped) we give it,
# we get differently named output files, which does not work well with snakemake.
# So we need to hardcode this here, it seems...
# We keep this function here, so that if we ever add an additinal step here, we will (hopefully)
# notice, and adapt this function as well.
def get_mapped_read_infix():
    # Same logic as the function above.
    result = "sorted"
    if config["settings"]["filter-mapped-reads"]:
        result = "filtered"
    if config["settings"]["clip-read-overlaps"]:
        result = "clipped"
    return result

# =================================================================================================
#     Mark Duplicates
# =================================================================================================

# Switch to the chosen duplicate marker tool
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

elif config["settings"]["duplicates-tool"] == "dedup":

    # Use `dedup`
    include: "duplicates-dedup.smk"

else:
    raise Exception("Unknown duplicates-tool: " + config["settings"]["duplicates-tool"])

# =================================================================================================
#     Base Quality Score Recalibration
# =================================================================================================

if config["settings"]["recalibrate-base-qualities"]:
    include: "mapping-recalibrate.smk"

# =================================================================================================
#     Indexing
# =================================================================================================

rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "logs/samtools/index/{prefix}.log"
    group:
        "mapping_extra"
    wrapper:
        "0.51.3/bio/samtools/index"

# =================================================================================================
#     Final Mapping Result
# =================================================================================================

# At this point, we have several choices of which files we want to hand down to the next
# pipleine step. We offer a function so that downstream does not have to deal with this.
# The way this works is as follows: Each optional step that could be applied
# overwrites the previous setting, until we have the final file that we want to use downstream.
# As each of these optional steps itself also takes care of using previous optional steps
# (see above), what we get here is the file name of the last step that we want to apply.

def get_mapping_result(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}-{unit}.sorted.bam"

    # case 2: filtering via samtools view
    if config["settings"]["filter-mapped-reads"]:
        f = "mapped/{sample}-{unit}.filtered.bam"

    # case 3: clipping reads with BamUtil
    if config["settings"]["clip-read-overlaps"]:
        f = "mapped/{sample}-{unit}.clipped.bam"

    # case 4: remove duplicates
    if config["settings"]["remove-duplicates"]:
        f = "dedup/{sample}-{unit}.bam"

    # case 5: recalibrate base qualities
    if config["settings"]["recalibrate-base-qualities"]:
        f = "recal/{sample}-{unit}.bam"

    # Additionally, this function is run for getting bai files as well
    if bai:
        f += ".bai"

    return f

# Return the bam file(s) for a given sample.
# Get all aligned reads of given sample, with all its units.
# This is where all units are merged together. The function also automatically gets
# which of the mapping results to use, depending on the config setting (whether to remove
# duplicates, and whether to recalibrate the base qualities), by using the get_mapping_result
# function, that gives the respective files depending on the config.
def get_sample_bams(sample):
    return expand(
        get_mapping_result(),
        sample=sample,
        unit=get_sample_units(sample)
        # unit=config["global"]["samples"].loc[sample].unit
    )

# Return the bai file(s) for a given sample
def get_sample_bais(sample):
    return expand(
        get_mapping_result(True),
        sample=sample,
        unit=get_sample_units(sample)
        # unit=config["global"]["samples"].loc[sample].unit
    )

# Return the bam file(s) for a sample, given a wildcard param from a rule.
def get_sample_bams_wildcard(wildcards):
    return get_sample_bams( wildcards.sample )

# Return the bai file(s) for a sample, given a wildcard param from a rule.
def get_sample_bais_wildcard(wildcards):
    return get_sample_bais( wildcards.sample )

# Return the bam file(s) for all samples
def get_all_bams():
    # Make a list of all bams in the order as the samples list.
    res = list()
    for su in config["global"]["sample-units"]:
        res.append( get_mapping_result().format( sample=su[0], unit=su[1] ))
    return res

    # The below approach gives the bams in sample-first order, which we do not want.
    # return list(chain.from_iterable( [
    #     get_sample_bams(sample) for sample in config["global"]["sample-names"]
    # ] ))

# Return the bai file(s) for all samples
def get_all_bais():
    res = list()
    for su in config["global"]["sample-units"]:
        res.append( get_mapping_result(True).format( sample=su[0], unit=su[1] ))
    return res

    # return list(chain.from_iterable( [
    #     get_sample_bais(sample) for sample in config["global"]["sample-names"]
    # ] ))

# Simple wildcard resolution.
def get_sample_bams_wildcards(wildcards):
    return get_sample_bams(wildcards.sample)

# =================================================================================================
#     All bams, but not SNP calling
# =================================================================================================

# This alternative target rule executes all steps up to th mapping, and yields the final bam
# files that would otherwise be used for variant calling in the downstream process.
# That is, depending on the config, these are the sorted, filtered, remove duplicates, or
# recalibrated base qualities bam files.
rule all_bams:
    input:
        get_all_bams()

# The `all_bams` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all_bams
