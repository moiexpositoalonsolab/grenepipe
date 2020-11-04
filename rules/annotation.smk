# =================================================================================================
#     SnpEff Setup
# =================================================================================================

# We store the snpeff database in the directory where the reference genome is.
# This way, we do not have to download it again and again for different runs.
def get_snpeff_db_path():
    return os.path.join(os.path.dirname( config["data"]["reference"]["genome"] ), "snpeff-db" ) + "/"

# We separate download from usage, so that we can better see progress and errors.
rule snpeff_db:
    output:
        # wildcard {reference} may be anything listed in the first column of `snpeff databases`
        directory(get_snpeff_db_path() + "{reference}")
    log:
        "logs/snpeff-download-{reference}.log"
    group:
        "snpeff"
    params:
        reference="{reference}"
    wrapper:
        "0.55.1/bio/snpeff/download"

# Rule is not submitted as a job to the cluster.
localrules: snpeff_db

# =================================================================================================
#     SnpEff
# =================================================================================================

rule snpeff:
    input:
        # (vcf, bcf, or vcf.gz)
        calls="filtered/all.vcf.gz",

        # path to reference db downloaded with the snpeff download wrapper above
        db=get_snpeff_db_path() + config["data"]["reference"]["name"]
    output:
        # annotated calls (vcf, bcf, or vcf.gz)
        calls=report("annotated/all.vcf.gz", caption="../reports/vcf.rst", category="Calls"),

        # summary statistics (in HTML), optional
        stats=report("snpeff/all.html", category="Calls"),

        # summary statistics in CSV, optional
        csvstats="snpeff/all.csv"
    log:
        "logs/snpeff.log"
    group:
        "snpeff"
    params:
        # optional parameters (e.g., max memory 4g)
        # For finding the chromosome names used by snpeff, add `-v` here
        extra=config["params"]["snpeff"]["extra"]
    wrapper:
        "0.55.1/bio/snpeff/annotate"

# Rule is not submitted as a job to the cluster.
localrules: snpeff
