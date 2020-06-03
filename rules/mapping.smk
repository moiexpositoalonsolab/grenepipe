from itertools import chain

# =================================================================================================
#     Read Mapping
# =================================================================================================

# Switch to the chosen mapper
if config["settings"]["mapping-tool"] == "bwamem":

    # Use `bwa mem`
    include: "mapping-bwa-mem.smk"

elif config["settings"]["mapping-tool"] == "bowtie2":

    # Use `bowtie2`
    include: "mapping-bowtie2.smk"

else:
    raise Exception("Unknown mapping-tool: " + config["settings"]["mapping-tool"])

# =================================================================================================
#     Mark Duplicates
# =================================================================================================

rule mark_duplicates:
    input:
        config["rundir"] + "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=config["rundir"] + "dedup/{sample}-{unit}.bam",
        metrics=config["rundir"] + "qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        config["rundir"] + "logs/picard/dedup/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/picard/dedup/{sample}-{unit}.bench.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    group:
        "mapping"
    wrapper:
        "0.51.3/bio/picard/markduplicates"

# =================================================================================================
#     Base Quality Score Recalibration
# =================================================================================================

def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = config["rundir"] + "mapped/{sample}-{unit}.sorted.bam"

    if config["settings"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = config["rundir"] + "dedup/{sample}-{unit}.bam"

    if bai:
        if config["settings"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
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
        known=config["data"]["reference"].get("known-variants") # empty if key not present
    output:
        bam=protected(config["rundir"] + "recal/{sample}-{unit}.bam")
    params:
        extra=get_gatk_regions_param() + " " + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        config["rundir"] + "logs/gatk/bqsr/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/gatk/bqsr/{sample}-{unit}.bench.log"
    group:
        "mapping"
    wrapper:
        "0.51.3/bio/gatk/baserecalibrator"

# =================================================================================================
#     Indexing
# =================================================================================================

rule bam_index:
    input:
        config["rundir"] + "{prefix}.bam"
    output:
        config["rundir"] + "{prefix}.bam.bai"
    log:
        config["rundir"] + "logs/samtools/index/{prefix}.log"
    group:
        "mapping"
    wrapper:
        "0.51.3/bio/samtools/index"

# =================================================================================================
#     Final Mapping Result
# =================================================================================================

# At this point, we have several choices of which files we want to hand down to the next
# pipleine step. We offer a function so that downstream does not have to deal with this.

def get_mapping_result(bai=False):
    # case 1: no duplicate removal
    f = config["rundir"] + "mapped/{sample}-{unit}.sorted.bam"

    # case 2: remove duplicates
    if config["settings"]["remove-duplicates"]:
        f = config["rundir"] + "dedup/{sample}-{unit}.bam"

    # case 3: recalibrate base qualities
    if config["settings"]["recalibrate-base-qualities"]:
        f = config["rundir"] + "recal/{sample}-{unit}.bam"

    # Additionally, this function is run for getting bai files as well
    if bai:
        f += ".bai"

    return f

# Return the bam file(s) for a given sample
def get_sample_bams(sample):
    return expand(get_mapping_result(), sample=sample, unit=samples.loc[sample].unit)

# Return the bai file(s) for a given sample
def get_sample_bais(sample):
    return expand(get_mapping_result(True), sample=sample, unit=samples.loc[sample].unit)

# Return the bam file(s) for all samples
def get_all_bams():
    return list(chain.from_iterable( [ get_sample_bams(sample) for sample in sample_names ] ))

# Return the bai file(s) for all samples
def get_all_bais():
    return list(chain.from_iterable( [ get_sample_bais(sample) for sample in sample_names ] ))
