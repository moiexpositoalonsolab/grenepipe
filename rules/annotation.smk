# =================================================================================================
#     SnpEff
# =================================================================================================

rule snpeff:
    input:
        config["rundir"] + "filtered/all.vcf.gz"
    output:
        vcf=report(config["rundir"] + "annotated/all.vcf.gz", caption="../reports/vcf.rst", category="Calls"),
        csvstats=config["rundir"] + "snpeff/all.csv"
    log:
        config["rundir"] + "logs/snpeff.log"
    params:
        reference=config["data"]["reference"]["name"],
        extra="-Xmx6g"
    wrapper:
        "0.27.1/bio/snpeff"
