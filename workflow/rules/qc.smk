# We use a switch to decide when we are running the pipeline normally with a set of fastq files,
# and when we are running it instead with a set of bam files to start with.
# In the latter case, some of the fastq and bam qc tools are not applicable.
use_fastq_input = not ("mappings-table" in config["data"] and config["data"]["mappings-table"])

# We outsource the individual QC steps for clarity.
if use_fastq_input:

    include: "qc-fastq.smk"


include: "qc-bam.smk"
include: "qc-vcf.smk"


# =================================================================================================
#     MultiQC
# =================================================================================================


# Unfortunately, in some cluster environments, multiqc does not work due to char encoding issues, see
# https://github.com/ewels/MultiQC/issues/484 ... If you run into this issue, try running it locally.
#
# We use the output touch files here only, in order to keep the snakemake log a bit cleaner.
# In the past, we used _all_ files of all QC tools directly as input here, which lead to the multiqc
# rule in the snakemake log be incredibly long. However, the multiqc wrapper that we are using uses
# the directories of these files anyway, so we can just start with these in the first place.
# That is, our `done` files are in the same dir as the actual QC files, so multiqc still finds them.


rule multiqc:
    input:
        # Fastq QC tools
        "qc/fastqc/fastqc.done" if use_fastq_input else [],
        "trimming/trimming-reports.done" if use_fastq_input else [],
        # Mapping QC tools
        "qc/samtools-stats/samtools-stats.done",
        "qc/samtools-flagstat/samtools-flagstat.done",
        "qc/qualimap/qualimap.done",
        "qc/picard/collectmultiplemetrics.done",
        get_dedup_done() if use_fastq_input and config["settings"]["remove-duplicates"] else [],
        # VCF QC tools, if requested
        # We only need the stats.vchk file of bcftools stats, but request the pdf here as well,
        # so that the bctools internal plots get generated by the bcftools_stats_plot rule as well.
        "qc/bcftools-stats/stats.vchk" if config["settings"]["bcftools-stats"] else [],
        "qc/bcftools-stats/summary.pdf" if config["settings"]["bcftools-stats"] else [],
        # Annotations, if requested
        "annotation/snpeff.csv" if config["settings"]["snpeff"] else [],
        "annotation/vep_summary.html" if config["settings"]["vep"] else [],
        # Damage profiling, if requested
        "damage/mapdamage/mapdamage.done"
        if use_fastq_input and config["settings"]["mapdamage"]
        else [],
        "damage/damageprofiler/damageprofiler.done"
        if use_fastq_input and config["settings"]["damageprofiler"]
        else [],
    output:
        report("qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control"),
        "qc/multiqc.zip",
    params:
        config["params"]["multiqc"]["extra"],
    log:
        "logs/qc/multiqc.log",
    benchmark:
        "benchmarks/qc/multiqc.log"
    conda:
        # We use a conda environment on top of the wrapper, as the wrapper always causes
        # issues with missing python modules and mismatching program versions and stuff...
        "../envs/multiqc.yaml"
    wrapper:
        "v3.13.6/bio/multiqc"


# script:
#     # We use our own version of the wrapper here, to troubleshoot dependecy issues...
#     "../scripts/multiqc.py"

# Rule is not submitted as a job to the cluster.
# Edit: It is now, as it turns out to be quite the process for large datasets...
# localrules: multiqc

# =================================================================================================
#     All QC, but not SNP calling
# =================================================================================================


# This alternative target rule executes all quality control (QC) steps of read trimming and mapping,
# but does not call SNPs, and does not call snpeff. The result is mainly the MultiQC report (without
# the snpeff part however), as well as the fastqc reports.
rule all_qc:
    input:
        # Quality control
        "qc/multiqc.html",
        "mapping/final.done",
        # Reference genome statistics, provided in the prep steps
        config["data"]["reference-genome"] + ".seqkit",


# The `all_qc` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules:
    all_qc,
