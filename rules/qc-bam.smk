import platform

# =================================================================================================
#     samtools_stats
# =================================================================================================

rule samtools_stats:
    input:
        get_mapping_result()
    output:
        "qc/samtools-stats/{sample}-{unit}.txt"
    log:
        "logs/samtools-stats/{sample}-{unit}.log"
    benchmark:
        "benchmarks/samtools-stats/{sample}-{unit}.bench.log"
    group:
        "qc"
    wrapper:
        "0.27.1/bio/samtools/stats"

rule samtools_stats_collect:
    input:
        expand(
            "qc/samtools-stats/{u.sample}-{u.unit}.txt",
            u=config["global"]["samples"].itertuples()
        )
    output:
        touch("qc/samtools-stats/samtools-stats.done")

localrules: samtools_stats_collect

# =================================================================================================
#     samtools_flagstat
# =================================================================================================

rule samtools_flagstat:
    input:
        get_mapping_result()
    output:
        "qc/samtools-flagstats/{sample}-{unit}.txt"
    log:
        "logs/samtools-flagstats/{sample}-{unit}.log"
    benchmark:
        "benchmarks/samtools-flagstats/{sample}-{unit}.bench.log"
    group:
        "qc"
    wrapper:
        "0.64.0/bio/samtools/flagstat"

rule samtools_flagstat_collect:
    input:
        expand(
            "qc/samtools-flagstats/{u.sample}-{u.unit}.txt",
            u=config["global"]["samples"].itertuples()
        )
    output:
        touch("qc/samtools-flagstats/samtools-flagstats.done")

localrules: samtools_flagstat_collect

# =================================================================================================
#     qualimap
# =================================================================================================

# We might want to run qualimap on either each bam file per sample-unit individually,
# or on a merged bam file for a whole sample, merging all its units.
# The merging rule here is duplicated from pileup.smk and from HAF-pipe in frequency.smk,
# but for now, we keep it that way, in order to be more flexible with the specifics.
rule qualimap_merge_unit_bams:
    input:
        get_sample_bams_wildcards # provided in mapping.smk
    output:
        # Need to pick a different directory than the future output directory
        # of the main qualimap rule here, so as not to confuse snakemake...
        temp("qc/qualimap-bams/{sample}.merged.bam"),
        touch("qc/qualimap-bams/{sample}.merged.done")
    params:
        config["params"]["samtools"]["merge"]
    threads:
        config["params"]["samtools"]["merge-threads"]
    log:
        "logs/samtools/qualimap/merge-{sample}.log"
    group:
        "qualimap"
    wrapper:
        "0.74.0/bio/samtools/merge"

# Rule for running qualimap on each bam file of a sample and its units separately.
# If things of how qualimap is run need to be adjusted here, they will probably also need to be
# adjusted below in the rule that runs per sample with merged units.
rule qualimap_sample_unit:
    input:
        get_mapping_result()
    output:
        # The output of this rule (for now) is a directory, so we have no way of telling whether
        # the rule succeeded by that. We hence add the `done` file here as well, which (according
        # to the snakemake docs) will only be created after the successful execution of the rule.
        # We had instances in the past where cluster jobs of this rule failed, and left a
        # half-finished directory behind; let's hope that this avoids that.
        outdir=directory("qc/qualimap/{sample}-{unit}"),
        done=touch("qc/qualimap/{sample}-{unit}.done")

        # Alternatively, specify all individual files, which gives more control,
        # but also more spammed output
        # "qc/qualimap/{sample}-{unit}/genome_results.txt",
        # "qc/qualimap/{sample}-{unit}/qualimapReport.html",
        # "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/coverage_histogram.txt",
        # "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/genome_fraction_coverage.txt",
        # "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
    params:
        extra=config["params"]["qualimap"]["extra"],
        outdir="qc/qualimap/{sample}-{unit}"
    threads:
        config["params"]["qualimap"]["threads"]
    log:
        "logs/qualimap/{sample}-{unit}_qualimap.log"
    group:
        "qualimap"
    conda:
        "../envs/qualimap.yaml"
    shell:
        "unset DISPLAY; qualimap bamqc -bam {input} -nt {threads} "
        "-outdir {params.outdir} -outformat HTML "
        "{params.extra} > {log} 2>&1"

# Rule for running qualimap for each bam file of a sample, merging its units.
# This is basically the same as above, but using the merged one instead.
# We are not repeating all comments of above here, to keep it simple. See above for details.
rule qualimap_sample_merged:
    input:
        "qc/qualimap-bams/{sample}.merged.bam"
    output:
        outdir=directory("qc/qualimap/{sample}"),
        done=touch("qc/qualimap/{sample}.done")
    params:
        extra=config["params"]["qualimap"]["extra"],
        outdir="qc/qualimap/{sample}"
    threads:
        config["params"]["qualimap"]["threads"]
    log:
        "logs/qualimap/{sample}_qualimap.log"
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
        # Get either the merged per-sample qualimap files, or the per-sample and per-unit ones.
        expand(
            "qc/qualimap/{sample}.done",
            sample=config["global"]["sample-names"]
        ) if config["params"]["qualimap"].get("merge-units", True) else
        expand(
            "qc/qualimap/{u.sample}-{u.unit}.done",
            u=config["global"]["samples"].itertuples()
        )

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
        touch("qc/qualimap/qualimap.done")

localrules: qualimap_collect

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

rule picard_collectmultiplemetrics:
    input:
        bam=get_mapping_result(),
        ref=config["data"]["reference-genome"]
    output:
        expand( "qc/picard/{{sample}}-{{unit}}{ext}", ext=picard_collectmultiplemetrics_exts())
    log:
        "logs/picard/multiple_metrics/{sample}-{unit}.log"
    params:
        config["params"]["picard"]["CollectMultipleMetrics"]["extra"] + (
            " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
            if platform.system() == "Darwin"
            else ""
        )
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
            "qc/picard/{u.sample}-{u.unit}{ext}",
            u=config["global"]["samples"].itertuples(),
            ext=picard_collectmultiplemetrics_exts()
        )
    output:
        touch("qc/picard/collectmultiplemetrics.done")

localrules: picard_collectmultiplemetrics_collect

# =================================================================================================
#     Deduplication
# =================================================================================================

# Different dedup tools produce different summary files. See above for details.
# This function simply returns these file names as strings, without replacing the wildcards.
# Then, when the function is called below, these are expanded.
def get_dedup_report():
    # Switch to the chosen duplicate marker tool
    if config["settings"]["duplicates-tool"] == "picard":
        return "qc/dedup/{u.sample}-{u.unit}.metrics.txt"
    elif config["settings"]["duplicates-tool"] == "dedup":
        return "dedup/{u.sample}-{u.unit}.dedup.json"
    else:
        raise Exception("Unknown duplicates-tool: " + config["settings"]["duplicates-tool"])

def get_dedup_done():
    # Switch to the chosen duplicate marker tool
    if config["settings"]["duplicates-tool"] == "picard":
        return "qc/dedup/picard.done"
    elif config["settings"]["duplicates-tool"] == "dedup":
        return "dedup/dedup.done"
    else:
        raise Exception("Unknown duplicates-tool: " + config["settings"]["duplicates-tool"])

rule dedup_reports_collect:
    input:
        expand(
            get_dedup_report(),
            u=config["global"]["samples"].itertuples()
        )
    output:
        touch( get_dedup_done() )

localrules: dedup_reports_collect
