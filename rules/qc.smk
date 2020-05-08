# =================================================================================================
#     Basic QC Stats
# =================================================================================================

rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html=config["rundir"] + "qc/fastqc/{sample}-{unit}.html",
        zip=config["rundir"] + "qc/fastqc/{sample}-{unit}.zip"
    log:
        config["rundir"] + "logs/fastqc/{sample}-{unit}.log"
    wrapper:
        "0.27.1/bio/fastqc"

rule samtools_stats:
    input:
        get_mapping_result()
    output:
        config["rundir"] + "qc/samtools-stats/{sample}-{unit}.txt"
    log:
        config["rundir"] + "logs/samtools-stats/{sample}-{unit}.log"
    wrapper:
        "0.27.1/bio/samtools/stats"

# =================================================================================================
#     MultiQC
# =================================================================================================

rule multiqc:
    input:
        expand(config["rundir"] + "qc/samtools-stats/{u.sample}-{u.unit}.txt", u=samples.itertuples()),
        expand(config["rundir"] + "qc/fastqc/{u.sample}-{u.unit}.zip", u=samples.itertuples()),
        expand(config["rundir"] + "qc/dedup/{u.sample}-{u.unit}.metrics.txt", u=samples.itertuples())

        # "snpeff/all.csv" TODO
    output:
        report(config["rundir"] + "qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        config["rundir"] + "logs/multiqc.log"
    wrapper:
        "0.55.1/bio/multiqc"
