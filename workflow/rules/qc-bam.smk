import platform


# The bam QC rules can either be run on the samples right after merging their units,
# or after the further downstream processing steps (filtering, clipping, dedup, recal etc)
# have been run. Here, we offer a simple function to switch between these inputs,
# in order to avoid duplicating thise code below.
# It expects the tool and key from the config.yaml
def bam_qc_input(tool, key, wildcards):
    if key not in config["params"][tool]:
        raise Exception("config key " + key + " not found for tool " + tool)
    if "mappings-table" in config["data"] and config["data"]["mappings-table"]:
        return get_sample_bams(wildcards.sample)
    elif config["params"][tool][key] == "processed":
        return get_sample_bams(wildcards.sample)
    elif config["params"][tool][key] == "merged":
        return "mapping/merged/{sample}.bam".format(sample=wildcards.sample)
    else:
        raise Exception(
            "Unknown setting for " + tool + " " + key + ": " + config["params"][tool][key]
        )


# =================================================================================================
#     samtools stats and flagstat
# =================================================================================================


rule samtools_stats:
    input:
        # bam_qc_input("samtools", "stats-bams"),
        lambda wildcards: bam_qc_input("samtools", "stats-bams", wildcards),
    output:
        "qc/samtools-stats/{sample}.txt",
    log:
        "logs/qc/samtools-stats/{sample}.log",
    benchmark:
        "benchmarks/qc/samtools-stats/{sample}.log"
    group:
        "qc"
    conda:
        "../envs/samtools.yaml"
    wrapper:
        "0.27.1/bio/samtools/stats"


rule samtools_stats_collect:
    input:
        expand("qc/samtools-stats/{sample}.txt", sample=config["global"]["sample-names"]),
    output:
        touch("qc/samtools-stats/samtools-stats.done"),


localrules:
    samtools_stats_collect,


rule samtools_flagstat:
    input:
        # bam_qc_input("samtools", "flagstat-bams"),
        lambda wildcards: bam_qc_input("samtools", "flagstat-bams", wildcards),
    output:
        "qc/samtools-flagstat/{sample}.txt",
    log:
        "logs/qc/samtools-flagstat/{sample}.log",
    benchmark:
        "benchmarks/qc/samtools-flagstat/{sample}.log"
    group:
        "qc"
    conda:
        "../envs/samtools.yaml"
    wrapper:
        "0.64.0/bio/samtools/flagstat"


rule samtools_flagstat_collect:
    input:
        expand("qc/samtools-flagstat/{sample}.txt", sample=config["global"]["sample-names"]),
    output:
        touch("qc/samtools-flagstat/samtools-flagstat.done"),


localrules:
    samtools_flagstat_collect,


# =================================================================================================
#     qualimap
# =================================================================================================


# Rule for running qualimap for each bam file of a sample, with its units merged.
rule qualimap_sample:
    input:
        # bam_qc_input("qualimap", "bams"),
        lambda wildcards: bam_qc_input("qualimap", "bams", wildcards),
    output:
        # The output of this rule (for now) is a directory, so we have no way of telling whether
        # the rule succeeded by that. We hence add the `done` file here as well, which (according
        # to the snakemake docs) will only be created after the successful execution of the rule.
        # We had instances in the past where cluster jobs of this rule failed, and left a
        # half-finished directory behind; let's hope that this avoids that.
        outdir=directory("qc/qualimap/{sample}"),
        done=touch("qc/qualimap/{sample}.done"),
        # Alternatively, specify all individual files, which gives more control,
        # but also more spammed output
        # "qc/qualimap/{sample}/genome_results.txt",
        # "qc/qualimap/{sample}/qualimapReport.html",
        # "qc/qualimap/{sample}/raw_data_qualimapReport/coverage_histogram.txt",
        # "qc/qualimap/{sample}/raw_data_qualimapReport/genome_fraction_coverage.txt",
        # "qc/qualimap/{sample}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
    params:
        extra=config["params"]["qualimap"]["extra"],
        outdir="qc/qualimap/{sample}",
    threads: config["params"]["qualimap"]["threads"]
    log:
        "logs/qc/qualimap/{sample}_qualimap.log",
    group:
        "qualimap"
    conda:
        "../envs/qualimap.yaml"
    shell:
        "unset DISPLAY; qualimap bamqc -bam {input} -nt {threads} "
        "-outdir {params.outdir} -outformat HTML "
        "{params.extra} > {log} 2>&1"


rule qualimap_collect:
    input:
        expand("qc/qualimap/{sample}.done", sample=config["global"]["sample-names"]),
        # Old, explicitly request each output file. Also still uses units...
        # expand(
        #     "qc/qualimap/{u.sample}-{u.unit}/genome_results.txt",
        #     u=config["global"]["samples"].itertuples()
        # ),
        # expand(
        #     "qc/qualimap/{u.sample}-{u.unit}/qualimapReport.html",
        #     u=config["global"]["samples"].itertuples()
        # ),
        # expand(
        #     "qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/coverage_histogram.txt",
        #     u=config["global"]["samples"].itertuples()
        # ),
        # expand(
        #     "qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/genome_fraction_coverage.txt",
        #     u=config["global"]["samples"].itertuples()
        # ),
        # expand(
        #     "qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
        #     u=config["global"]["samples"].itertuples()
        # ),
    output:
        touch("qc/qualimap/qualimap.done"),


localrules:
    qualimap_collect,


# =================================================================================================
#     Picard CollectMultipleMetrics
# =================================================================================================


# The snakemake wrapper for picard/collectmultiplemetrics uses output file extensions
# to select the different tools for the metrics.
# Usable extensions (and which tools they implicitly call) are listed here:
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/collectmultiplemetrics.html
def picard_collectmultiplemetrics_exts():
    res = []
    if config["params"]["picard"]["CollectMultipleMetrics"]["AlignmentSummaryMetrics"]:
        res.append(".alignment_summary_metrics")
    if config["params"]["picard"]["CollectMultipleMetrics"]["BaseDistributionByCycle"]:
        res.append(".base_distribution_by_cycle_metrics")
        res.append(".base_distribution_by_cycle.pdf")
    if config["params"]["picard"]["CollectMultipleMetrics"]["GcBiasMetrics"]:
        res.append(".gc_bias.detail_metrics")
        res.append(".gc_bias.summary_metrics")
        res.append(".gc_bias.pdf")
    if config["params"]["picard"]["CollectMultipleMetrics"]["InsertSizeMetrics"]:
        res.append(".insert_size_metrics")
        res.append(".insert_size_histogram.pdf")
    if config["params"]["picard"]["CollectMultipleMetrics"]["QualityByCycleMetrics"]:
        res.append(".quality_by_cycle_metrics")
        res.append(".quality_by_cycle.pdf")
    if config["params"]["picard"]["CollectMultipleMetrics"]["QualityScoreDistributionMetrics"]:
        res.append(".quality_distribution_metrics")
        res.append(".quality_distribution.pdf")
    if config["params"]["picard"]["CollectMultipleMetrics"]["QualityYieldMetrics"]:
        res.append(".quality_yield_metrics")
    # if config["params"]["picard"]["CollectMultipleMetrics"]["RnaSeqMetrics"]:
    #     res.append(".rna_metrics")
    return res


# Need a function here, to handle wildcards, as a lambda with a named input does not work...
# def picard_collectmultiplemetrics_input(wildcards):
#     return bam_qc_input("picard", "CollectMultipleMetrics-bams", wildcards),


rule picard_collectmultiplemetrics:
    input:
        # bam=bam_qc_input("picard", "CollectMultipleMetrics-bams"),
        bam=lambda wildcards: bam_qc_input("picard", "CollectMultipleMetrics-bams", wildcards),
        # bam = picard_collectmultiplemetrics_input
        ref=config["data"]["reference-genome"],
    output:
        expand("qc/picard/{{sample}}{ext}", ext=picard_collectmultiplemetrics_exts()),
    log:
        "logs/qc/picard-collectmultiplemetrics/{sample}.log",
    params:
        config["params"]["picard"]["CollectMultipleMetrics-java-opts"]
        + " "
        + config["params"]["picard"]["CollectMultipleMetrics-extra"]
        + (" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true" if platform.system() == "Darwin" else ""),
    conda:
        "../envs/picard.yaml"
    script:
        # We use our own version of the wrapper here, which fixes issues with missing files in cases
        # where Picard does not have enough data for a specific metric to run.
        "../scripts/picard-collectmultiplemetrics.py"


# wrapper:
#     "0.72.0/bio/picard/collectmultiplemetrics"


rule picard_collectmultiplemetrics_collect:
    input:
        expand(
            "qc/picard/{sample}{ext}",
            sample=config["global"]["sample-names"],
            ext=picard_collectmultiplemetrics_exts(),
        ),
    output:
        touch("qc/picard/collectmultiplemetrics.done"),


localrules:
    picard_collectmultiplemetrics_collect,


# =================================================================================================
#     Deduplication
# =================================================================================================


# Different dedup tools produce different summary files. See above for details.
# This function simply returns these file names as strings, without replacing the wildcards.
# Then, when the function is called below, these are expanded.
def get_dedup_report():
    # Switch to the chosen duplicate marker tool
    if config["settings"]["duplicates-tool"] == "picard":
        return "qc/dedup/{sample}.metrics.txt"
    elif config["settings"]["duplicates-tool"] == "dedup":
        return "mapping/dedup/{sample}.dedup.json"
    else:
        raise Exception("Unknown duplicates-tool: " + config["settings"]["duplicates-tool"])


def get_dedup_done():
    # Switch to the chosen duplicate marker tool
    if config["settings"]["duplicates-tool"] == "picard":
        return "qc/dedup/picard.done"
    elif config["settings"]["duplicates-tool"] == "dedup":
        return "mapping/dedup/dedup.done"
    else:
        raise Exception("Unknown duplicates-tool: " + config["settings"]["duplicates-tool"])


rule dedup_reports_collect:
    input:
        expand(get_dedup_report(), sample=config["global"]["sample-names"]),
    output:
        touch(get_dedup_done()),


localrules:
    dedup_reports_collect,
