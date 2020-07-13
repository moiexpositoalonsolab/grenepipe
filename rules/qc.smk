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

# Unfortunately, in some environments, multiqc does not work due to encoding issues, see
# https://github.com/ewels/MultiQC/issues/484
# Hence, we here do not use the wrapper, but instead call the command manually.
rule multiqc:
    input:
        expand("qc/samtools-stats/{u.sample}-{u.unit}.txt", u=samples.itertuples()),
        expand("qc/fastqc/{u.sample}-{u.unit}.zip", u=samples.itertuples()),
        expand("qc/dedup/{u.sample}-{u.unit}.metrics.txt", u=samples.itertuples()),
        "snpeff/all.csv"
    output:
        report("qc/multiqc.html", caption="../reports/multiqc.rst", category="Quality control")
    log:
        "logs/multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    wrapper:
        "0.55.1/bio/multiqc"

# Rule is not submitted as a job to the cluster.
localrules: multiqc
