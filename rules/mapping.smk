# =================================================================================================
#     Trimming
# =================================================================================================

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp(config["rundir"] + "trimmed/{sample}-{unit}.fastq.gz")
    params:
        extra="",
        **config["params"]["trimmomatic"]["se"]
    log:
        config["rundir"] + "logs/trimmomatic/{sample}-{unit}.log"
    wrapper:
        "0.51.3/bio/trimmomatic/se"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp(config["rundir"] + "trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp(config["rundir"] + "trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp(config["rundir"] + "trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp(config["rundir"] + "trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog=config["rundir"] + "trimmed/{sample}-{unit}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    log:
        config["rundir"] + "logs/trimmomatic/{sample}-{unit}.log"
    wrapper:
        "0.51.3/bio/trimmomatic/pe"

# =================================================================================================
#     Read Mapping
# =================================================================================================

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if is_single_end(**wildcards):
        # single end sample
        return config["rundir"] + "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)
    else:
        # paired-end sample
        return expand(config["rundir"] + "trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards)

rule map_reads:
    input:
        reads=get_trimmed_reads
    output:
        config["rundir"] + "mapped/{sample}-{unit}.sorted.bam"
    log:
        config["rundir"] + "logs/bwa_mem/{sample}-{unit}.log"
    params:
        index=config["data"]["reference"]["genome"],

        # We need the read group tags, including `ID` and `SM`, as downstream tools use these.
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        # TODO add LD field as well for the unit?! http://www.htslib.org/workflow/

        sort="samtools",         # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate", # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 8
    wrapper:
        "0.51.3/bio/bwa/mem"

# TODO add threads: http://bio-bwa.sourceforge.net/bwa.shtml
# TODO add bowtie alternative

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
    params:
        config["params"]["picard"]["MarkDuplicates"]
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
