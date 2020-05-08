# =================================================================================================
#     Calls Table
# =================================================================================================

rule vcf_to_tsv:
    input:
        config["rundir"] + "annotated/all.vcf.gz"
    output:
        report(config["rundir"] + "tables/calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"

# =================================================================================================
#     Depth Stats Plot
# =================================================================================================

rule plot_stats:
    input:
        config["rundir"] + "tables/calls.tsv.gz"
    output:
        depths=report(config["rundir"] + "plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report(config["rundir"] + "plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot-depths.py"
