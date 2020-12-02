# =================================================================================================
#     Basic QC Stats
# =================================================================================================

def get_fastqc_input(wildcards):
    if config["params"]["fastqc"]["input"] == "samples":
        return get_fastq(wildcards)
    elif config["params"]["fastqc"]["input"] == "trimmed":
        return get_trimmed_reads(wildcards)
    else:
        raise Exception("Unknown fastqc input setting: " + config["params"]["fastqc"]["input"])

rule fastqc:
    input:
        unpack(get_fastqc_input)
    output:
        html="qc/fastqc/{sample}-{unit}_fastqc.html",
        zip="qc/fastqc/{sample}-{unit}_fastqc.zip"
    params:
        config["params"]["fastqc"]["extra"]
    log:
        "logs/fastqc/{sample}-{unit}.log"
    benchmark:
        "benchmarks/fastqc/{sample}-{unit}.bench.log"
    group:
        "qc"
    conda:
        "../envs/fastqc.yaml"
    script:
        # We use our own version of the wrapper here, as that wrapper is just badly implemented.
        "../scripts/fastqc.py"
    # wrapper:
    #     "0.27.1/bio/fastqc"

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

        # Alternatively, specify all individual files, which gives more control, but also more spammed ouutput
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
#     MultiQC
# =================================================================================================

# Different trimming tools produce different summary files. We here expand ourselves,
# because we need to retrieve the correct file type for each of them (single or paired end),
# which cannot easily be done with the simple snakemake expand function.
def get_trimming_reports():
    result=[]
    for smp in samples.itertuples():
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

# Unfortunately, in some environments, multiqc does not work due to char encoding issues, see
# https://github.com/ewels/MultiQC/issues/484 ... If you run into this issue, try running it locally.
rule multiqc:
    input:
        # Trimming
        get_trimming_reports(),

        # Dedup
        expand(get_dedup_report(), u=samples.itertuples()),

        # QC tools
        expand("qc/fastqc/{u.sample}-{u.unit}_fastqc.zip", u=samples.itertuples()),
        expand("qc/samtools-stats/{u.sample}-{u.unit}.txt", u=samples.itertuples()),
        expand("qc/samtools-flagstats/{u.sample}-{u.unit}.txt", u=samples.itertuples()),

        # Qualimap
        expand("qc/qualimap/{u.sample}-{u.unit}", u=samples.itertuples()),
        # expand("qc/qualimap/{u.sample}-{u.unit}/genome_results.txt", u=samples.itertuples()),
        # expand("qc/qualimap/{u.sample}-{u.unit}/qualimapReport.html", u=samples.itertuples()),
        # expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/coverage_histogram.txt", u=samples.itertuples()),
        # expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/genome_fraction_coverage.txt", u=samples.itertuples()),
        # expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt", u=samples.itertuples()),

        # Annotation
        "snpeff/all.csv" if config["settings"]["snpeff"] else [],

        # Damage
        expand("mapdamage/{u.sample}-{u.unit}/Runtime_log.txt", u=samples.itertuples())  if config["settings"]["mapdamage"] else [],
        expand("damageprofiler/{u.sample}-{u.unit}/DamageProfiler.log", u=samples.itertuples())  if config["settings"]["damageprofiler"] else []
    output:
        report("qc/multiqc.html", caption="../reports/multiqc.rst", category="Quality control")
    params:
        config["params"]["multiqc"]["extra"],
    log:
        "logs/multiqc.log"
    wrapper:
        "0.64.0/bio/multiqc"
    # conda:
    #     "../envs/multiqc.yaml"
    # script:
    #     # We use our own version of the wrapper here, to troubleshoot dependecy issues...
    #     "../scripts/multiqc.py"

# Rule is not submitted as a job to the cluster.
localrules: multiqc
