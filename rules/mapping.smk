from itertools import chain

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
#     Filtering Mapped Reads
# =================================================================================================

rule filter_mapped_reads:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        "mapped/{sample}-{unit}.filtered.bam"
    params:
        config["params"]["samtools"]["view"] + " -b"
    wrapper:
        "0.72.0/bio/samtools/view"

def get_mapped_reads(wildcards):
    """Get mapped reads of given sample-unit, either filtered or unfiltered"""
    if config["settings"]["filter-mapped-reads"]:
        return "mapped/{sample}-{unit}.filtered.bam".format(**wildcards)
    else:
        return "mapped/{sample}-{unit}.sorted.bam".format(**wildcards)

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

def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}-{unit}.sorted.bam"

    # case 2: filtering via samtools view
    if config["settings"]["filter-mapped-reads"]:
        f = "mapped/{sample}-{unit}.filtered.bam"

    if config["settings"]["remove-duplicates"]:
        # case 3: remove duplicates
        f = "dedup/{sample}-{unit}.bam"

    if bai:
        if config["settings"].get("restrict-regions"):
            # case 4: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 5: no index needed
            return []
    else:
        return f

def get_gatk_regions_param(regions=config["settings"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["settings"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default

rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),
        known=config["data"]["reference"]["known-variants"]
    output:
        bam=protected("recal/{sample}-{unit}.bam")
    params:
        extra=get_gatk_regions_param() + " " + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log"
    benchmark:
        "benchmarks/gatk/bqsr/{sample}-{unit}.bench.log"
    group:
        "mapping_extra"
    wrapper:
        "0.51.3/bio/gatk/baserecalibrator"

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

def get_mapping_result(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}-{unit}.sorted.bam"

    # case 2: filtering via samtools view
    if config["settings"]["filter-mapped-reads"]:
        f = "mapped/{sample}-{unit}.filtered.bam"

    # case 3: remove duplicates
    if config["settings"]["remove-duplicates"]:
        f = "dedup/{sample}-{unit}.bam"

    # case 4: recalibrate base qualities
    if config["settings"]["recalibrate-base-qualities"]:
        if not config["data"]["reference"].get("known-variants"):
            raise Exception(
                "Setting recalibrate-base-qualities can only be activated if a known-variants "
                "file is also provided. As far as we are aware, GATK BaseRecalibrator needs this."
            )
        f = "recal/{sample}-{unit}.bam"

    # Additionally, this function is run for getting bai files as well
    if bai:
        f += ".bai"

    return f

# Return the bam file(s) for a given sample
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
