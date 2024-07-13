# =================================================================================================
#     bcftools stats
# =================================================================================================


rule bcftools_stats:
    input:
        # we use the filtered file if a filtering is done, or the unfiltered if not.
        calls=(
            "calling/filtered-all.vcf.gz"
            if not config["settings"]["filter-variants"] == "none"
            else "calling/genotyped-all.vcf.gz"
        ),
    output:
        "qc/bcftools-stats/stats.vchk",
    log:
        "logs/qc/bcftools-stats.log",
    params:
        config["params"]["bcftools"]["stats"],
    conda:
        "../envs/bcftools.yaml"
    group:
        "bcftools-stats"
    wrapper:
        "v1.7.0/bio/bcftools/stats"


rule bcftools_stats_plot:
    input:
        "qc/bcftools-stats/stats.vchk",
    output:
        # We only request the PDF here, but other files are generated as well.
        "qc/bcftools-stats/summary.pdf",
    log:
        "logs/qc/bcftools-stats-plot.log",
    params:
        outdir="qc/bcftools-stats",
        extra=config["params"]["bcftools"]["stats-plot"],
    conda:
        "../envs/bcftools.yaml"
    group:
        "bcftools-stats"
    shell:
        # According to this thread: https://stackoverflow.com/a/69671413/4184258
        # there are issues with pdflatex from conda (because, why would software ever work?!).
        # So we work our way around this by working with tectonic instead as well,
        # in case that pdflatex fails. For this to work, we also need to stop the bash pipefail
        # (strict mode) that is activated by snakemake from interfering here, so that the failing
        # plot command does not stop the whole exection,
        # see https://stackoverflow.com/a/11231970/4184258 for that.
        # Lastly, in a freak way, somehow using a simple `cd` here seems to fail on our cluster...
        # So instead we specify the path to the file directly. So weird.
        "( plot-vcfstats --prefix {params.outdir} {params.extra} {input} &> {log} || true ) ; "
        "if [ -f {output} ]; then "
        '    echo "Success with pdflatex" >> {log} ; '
        "else"
        '    echo "Failed with pdflatex" >> {log} ; '
        '    echo "Trying tectonic instead" >> {log} ; '
        "    tectonic {params.outdir}/summary.tex >> {log} 2>&1 ; "
        "fi"
