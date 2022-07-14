# =================================================================================================
#     FastQC
# =================================================================================================

# Make a list of all files that we want to run FastQC on.
# We store them in the global config, so that the rule can access them.
# There are several combinations of cases of which files we want to run fastqc on.
# We can have single or paired end raw read sample data; we can have the trimmed reads, which
# can also be either single or paired end, but also merged pairs. We here define a function
# that takes care of all theses cases and returns a list of all the ones we want to produce.
# The MultiQC rule then uses this list to request all those files, and the below get_fastqc_input()
# function resolves them into the input file paths.

assert "fastqc" not in config["global"]
config["global"]["fastqc"] = pd.DataFrame( columns = [ "sample", "unit", "id", "file" ])
def add_fastqc_file( sample, unit, id, file ):
    config["global"]["fastqc"] = pd.concat([
        config["global"]["fastqc"],
        pd.DataFrame({
            "sample": [sample],
            "unit":   [unit],
            "id":     [id],
            "file":   [file]
        })
    ], ignore_index=True)

    # Deprecated way of using append instead of concat. We use concat now to avoid warnings
    # with newer pandas versions, but keep the below for reference.
    # config["global"]["fastqc"] = config["global"]["fastqc"].append({
    #     "sample": sample,
    #     "unit":   unit,
    #     "id":     id,
    #     "file":   file
    # }, ignore_index=True)

if config["params"]["fastqc"]["input"] == "samples":
    # Simple case: raw fastq files from the samples.
    for smp in config["global"]["samples"].itertuples():
        add_fastqc_file( smp.sample, smp.unit, "R1", smp.fq1 )
        if isinstance(smp.fq2, str):
            add_fastqc_file( smp.sample, smp.unit, "R2", smp.fq2 )
elif config["params"]["fastqc"]["input"] == "trimmed":
    # Trimmed files, which can come in more varieties.
    for smp in config["global"]["samples"].itertuples():
        # We use a fake wildcard to get the function call to work.
        wc = snakemake.io.Wildcards()
        wc.sample = smp.sample
        wc.unit = smp.unit
        trimmed = get_trimmed_reads(wc)

        # Now let's see if we have merged them or not, and add to our result accordingly.
        if config["settings"]["merge-paired-end-reads"]:
            assert len(trimmed) == 1
            add_fastqc_file( smp.sample, smp.unit, "trimmed-merged", trimmed[0] )
        else:
            # If not merged, it's either single end or paired end.
            assert len(trimmed) == 1 or len(trimmed) == 2
            add_fastqc_file( smp.sample, smp.unit, "trimmed-R1", trimmed[0] )
            if len(trimmed) == 2:
                add_fastqc_file( smp.sample, smp.unit, "trimmed-R2", trimmed[1] )
else:
    raise Exception("Unknown fastqc input setting: " + config["params"]["fastqc"]["input"])

def get_fastqc_input(wildcards):
    return config["global"]["fastqc"].loc[
        (config["global"]["fastqc"]['sample'] == wildcards.sample) &
        (config["global"]["fastqc"]['unit']   == wildcards.unit  ) &
        (config["global"]["fastqc"]['id']     == wildcards.id    ),
        ["file"]
    ].file

rule fastqc:
    input:
        get_fastqc_input
    output:
        html="qc/fastqc/{sample}-{unit}-{id}_fastqc.html",
        zip="qc/fastqc/{sample}-{unit}-{id}_fastqc.zip"
    params:
        config["params"]["fastqc"]["extra"]
    log:
        "logs/fastqc/{sample}-{unit}-{id}.log"
    benchmark:
        "benchmarks/fastqc/{sample}-{unit}-{id}.bench.log"
    group:
        "qc"
    conda:
        "../envs/fastqc.yaml"
    script:
        # We use our own version of the wrapper here, as that wrapper is just badly implemented.
        "../scripts/fastqc.py"
    # wrapper:
    #     "0.27.1/bio/fastqc"

# =================================================================================================
#     Basic Mapping QC Stats
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

rule qualimap:
    input:
        get_mapping_result()
    output:
        outdir=directory("qc/qualimap/{sample}-{unit}")

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
        stderr="logs/qualimap/{sample}-{unit}_qualimap.stderr",
        stdout="logs/qualimap/{sample}-{unit}_qualimap.stdout"
    group:
        "qc"
    conda:
        "../envs/qualimap.yaml"
    shell:
        "unset DISPLAY; qualimap bamqc -bam {input} -nt {threads} "
        "-outdir {params.outdir} -outformat HTML "
        "{params.extra} > {log.stdout} 2> {log.stderr}"

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
        ref=config["data"]["reference"]["genome"]
    output:
        expand( "qc/picard/{{sample}}-{{unit}}{ext}", ext=picard_collectmultiplemetrics_exts())
    log:
        "logs/picard/multiple_metrics/{sample}-{unit}.log"
    params:
        config["params"]["picard"]["CollectMultipleMetrics"]["extra"]
    conda:
        "../envs/picard.yaml"
    script:
        # We use our own version of the wrapper here, which fixes issues with missing files in cases
        # where Picard does not have enough data for a specific metric to run.
        "../scripts/picard-collectmultiplemetrics.py"
    # wrapper:
    #     "0.72.0/bio/picard/collectmultiplemetrics"

# =================================================================================================
#     bcftools stats
# =================================================================================================

rule bcftools_stats:
    input:
        "filtered/all.vcf.gz"
    output:
        "qc/bcftools-stats/stats.vchk"
    log:
        "logs/bcftools-stats/bcftools.stats.log"
    params:
        config["params"]["bcftools"]["stats"]
    conda:
        "../envs/bcftools.yaml"
    wrapper:
        "v1.7.0/bio/bcftools/stats"

rule bcftools_stats_plot:
    input:
        "qc/bcftools-stats/stats.vchk"
    output:
        # We only request the PDF here, but other files are generated as well.
        "qc/bcftools-stats/summary.pdf"
    log:
        "logs/bcftools-stats/bcftools.stats.plot.log"
    params:
        outdir="qc/bcftools-stats",
        extra=config["params"]["bcftools"]["stats-plot"]
    conda:
        "../envs/bcftools.yaml"
    shell:
        # According to this thread: https://stackoverflow.com/a/69671413/4184258
        # there are issues with pdflatex and conda (because, why would software ever work?!).
        # So we work our way around this by working with tectonic instead as well...
        # We also need to stop the bash pipefail (strict mode) that is activated by snakemake
        # from interfering here, so that the failing plot command does not stop the whole exection,
        # see https://stackoverflow.com/a/11231970/4184258 for that.
        "( plot-vcfstats --prefix {params.outdir} {params.extra} {input} &> {log} || true ) ; "
        "if [ ! -f \"{output}\" ]; then "
        "    cd {params.outdir} ; "
        "    tectonic summary.tex >> {log} 2>&1 ; "
        "fi"

# =================================================================================================
#     MultiQC
# =================================================================================================

# Different trimming tools produce different summary files. We here expand ourselves,
# because we need to retrieve the correct file type for each of them (single or paired end),
# which cannot easily be done with the simple snakemake expand function.
def get_trimming_reports():
    result=[]
    for smp in config["global"]["samples"].itertuples():
        # The get_trimming_report() function is part of each trimming tool rule file.
        # We here hence call the respective correct function for each tool.
        result.append( get_trimming_report( smp.sample, smp.unit ))

        # Now append the file for the sample to the result list
        # if config["settings"]["trimming-tool"] == "adapterremoval":
            # result.append( "trimmed/" + smp.sample + "-" + smp.unit + "-" + suffix + ".settings" )
        # elif config["settings"]["trimming-tool"] == "cutadapt":
            # result.append( "trimmed/" + smp.sample + "-" + smp.unit + ".qc-" + suffix + ".txt" )
        # elif config["settings"]["trimming-tool"] == "fastp":
        #     result.append( "trimmed/" + smp.sample + "-" + smp.unit + "-" + suffix + "-fastp.json" )
        # elif config["settings"]["trimming-tool"] == "skewer":
        #     result.append( "trimmed/" + smp.sample + "-" + smp.unit + "-" + suffix + "-trimmed.log" )
        # elif config["settings"]["trimming-tool"] == "trimmomatic":
        #     result.append( "trimmed/" + smp.sample + "-" + smp.unit + ".trimlog.log" )
        # else:
        #     raise Exception("Unknown trimming-tool: " + config["settings"]["trimming-tool"])
    return result

# Different dedup tools produce different summary files. See above for details.
# This function simply returns these file names as strings, without replacing the wildcards.
# Then, when the function is called below, these are expanded.
def get_dedup_report():
    # Switch to the chosen duplicate marker tool
    if config["settings"]["duplicates-tool"] == "picard":
        return "qc/dedup/{u.sample}-{u.unit}.metrics.txt"
    elif config["settings"]["duplicates-tool"] == "dedup":
        return "dedup/{u.sample}-{u.unit}.sorted.dedup.json"
    else:
        raise Exception("Unknown duplicates-tool: " + config["settings"]["duplicates-tool"])

# Unfortunately, in some cluster environments, multiqc does not work due to char encoding issues, see
# https://github.com/ewels/MultiQC/issues/484 ... If you run into this issue, try running it locally.
rule multiqc:
    input:
        # FastQC
        expand(
            "qc/fastqc/{u.sample}-{u.unit}-{u.id}_fastqc.html",
            u=config["global"]["fastqc"].itertuples()
        ),
        expand(
            "qc/fastqc/{u.sample}-{u.unit}-{u.id}_fastqc.zip",
            u=config["global"]["fastqc"].itertuples()
        ),

        # Trimming
        get_trimming_reports(),

        # Dedup
        expand(
            get_dedup_report(),
            u=config["global"]["samples"].itertuples()
        ),

        # Mapping QC tools
        expand(
            "qc/samtools-stats/{u.sample}-{u.unit}.txt",
            u=config["global"]["samples"].itertuples()
        ),
        expand(
            "qc/samtools-flagstats/{u.sample}-{u.unit}.txt",
            u=config["global"]["samples"].itertuples()
        ),

        # Qualimap
        expand(
            "qc/qualimap/{u.sample}-{u.unit}",
            u=config["global"]["samples"].itertuples()
        ),
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

        # Picard CollectMultipleMetrics
        expand(
            "qc/picard/{u.sample}-{u.unit}{ext}",
            u=config["global"]["samples"].itertuples(),
            ext=picard_collectmultiplemetrics_exts()
        ),

        # bcftools stats
        # We only require the stats.vchk file, but request the pdf here as well,
        # so that the bctools internal plots get generated by the bcftools_stats_plot rule as well.
        "qc/bcftools-stats/stats.vchk",
        "qc/bcftools-stats/summary.pdf",

        # Annotation
        "annotated/snpeff.csv" if config["settings"]["snpeff"] else [],
        "annotated/vep_summary.html" if config["settings"]["vep"] else [],

        # Damage
        expand(
            "mapdamage/{u.sample}-{u.unit}/Runtime_log.txt",
            u=config["global"]["samples"].itertuples()
        ) if config["settings"]["mapdamage"] else [],
        expand(
            "damageprofiler/{u.sample}-{u.unit}/DamageProfiler.log",
            u=config["global"]["samples"].itertuples()
        ) if config["settings"]["damageprofiler"] else []
    output:
        report("qc/multiqc.html", caption="../reports/multiqc.rst", category="Quality control")
    params:
        config["params"]["multiqc"]["extra"],
    log:
        "logs/multiqc.log"
    conda:
        # We use a conda environment on top of the wrapper, as the wrapper always causes
        # issues with missing python modules and mismatching program versions and stuff...
        "../envs/multiqc.yaml"
    wrapper:
        "0.74.0/bio/multiqc"
    # script:
    #     # We use our own version of the wrapper here, to troubleshoot dependecy issues...
    #     "../scripts/multiqc.py"

# Rule is not submitted as a job to the cluster.
# Edit: It is now, as it turns out to be quite the process for large datasets...
# localrules: multiqc
