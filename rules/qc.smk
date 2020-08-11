# =================================================================================================
#     Basic QC Stats
# =================================================================================================

rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}-{unit}.html",
        zip="qc/fastqc/{sample}-{unit}.zip"
    log:
        "logs/fastqc/{sample}-{unit}.log"
    benchmark:
        "benchmarks/fastqc/{sample}-{unit}.bench.log"
    group:
        "qc"
    wrapper:
        "0.27.1/bio/fastqc"

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

# =================================================================================================
#     MultiQC
# =================================================================================================

# Different dedup tools produce different summary files. This function simply returns these file
# names as strings, without replacing the wildcards. Then, when the function is called below,
# these are expanded.
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
        expand("qc/samtools-stats/{u.sample}-{u.unit}.txt", u=samples.itertuples()),
        expand("qc/fastqc/{u.sample}-{u.unit}.zip", u=samples.itertuples()),
        expand(get_dedup_report(), u=samples.itertuples()),
        "snpeff/all.csv"
    output:
        report("qc/multiqc.html", caption="../reports/multiqc.rst", category="Quality control")
    params:
        config["params"]["multiqc"]["extra"],
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    wrapper:
        "0.64.0/bio/multiqc"

# Rule is not submitted as a job to the cluster.
localrules: multiqc
