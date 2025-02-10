# =================================================================================================
#     Calls Table
# =================================================================================================


rule vcf_to_tsv:
    input:
        "annotation/snpeff.vcf.gz",
    output:
        report("tables/calls.tsv.gz", caption="../report/calls.rst", category="Calls"),
    log:
        "logs/stats/vcf_to_tsv.log",
    conda:
        "../envs/rbt.yaml"
    group:
        "stats"
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"


# TODO this rule fails for freebayes for some weird reason, see
# https://github.com/rust-bio/rust-bio-tools/issues/52#issuecomment-626289782
# fix it, once that github issue has an answer

# =================================================================================================
#     Depth Stats Plot
# =================================================================================================


rule plot_stats:
    input:
        "tables/calls.tsv.gz",
    output:
        depths=report("plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report("plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots"),
    log:
        "logs/stats/plot-depths.log",
    conda:
        "../envs/stats.yaml"
    group:
        "stats"
    script:
        "../scripts/plot-depths.py"


# =================================================================================================
#     Sequences per Sample
# =================================================================================================


rule seqs_per_sample:
    output:
        "tables/sample-sizes.tsv",
    params:
        samples=config["data"]["samples-table"],
    script:
        "../scripts/sample-sizes.py"


# =================================================================================================
#     Allele Frequency Table
# =================================================================================================


rule frequency_table:
    input:
        # we use the filtered file if a filtering is done, or the unfiltered if not.
        calls=(
            "calling/filtered-all.vcf.gz"
            if not config["settings"]["filter-variants"] == "none"
            else "calling/genotyped-all.vcf.gz"
        ),
        done=(
            "calling/filtered-all.vcf.gz.done"
            if not config["settings"]["filter-variants"] == "none"
            else "calling/genotyped-all.vcf.gz.done"
        ),
    output:
        "tables/frequencies.tsv",
    params:
        fields=config["settings"]["frequency-table-fields"],
    log:
        "logs/stats/frequency-table.log",
    conda:
        "../envs/frequency-table.yaml"
    script:
        "../scripts/frequency-table.py"
