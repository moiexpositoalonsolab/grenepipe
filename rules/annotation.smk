# =================================================================================================
#     SnpEff Setup
# =================================================================================================

# Get the path where we download the database to, so that we do not have to download it all the time.
# If not provided by the user, we use the directory where the reference genome is.
# Not used with a custom db, in which case we just use the provided dir directly.
def get_snpeff_db_download_path():
    if config["params"]["snpeff"]["download-dir"]:
        # Return path from config, and make sure that it has trailing slash
        return os.path.join( config["params"]["snpeff"]["download-dir"], '' )
    else:
        # Use the ref genome path, with trailing slash
        return os.path.join(
            os.path.dirname( config["data"]["reference-genome"] ), "snpeff-db"
        ) + "/"

# We separate download from usage, so that we can better see progress and errors,
# and to avoid unneccessary re-downloads.
rule snpeff_db:
    output:
        # wildcard {reference} may be anything listed in the first column of `snpeff databases`
        directory(get_snpeff_db_download_path() + "{reference}")
    log:
        "logs/snpeff-download-{reference}.log"
    group:
        "snpeff"
    params:
        reference="{reference}"
    conda:
        # As always, the wrapper does not specify that python, pandas, and numpy are needed,
        # and if those are misspecified on the cluster (as in our case), it just fails...
        # So let's provide a working combination of versions of these tools.
        "../envs/snpeff.yaml"
    wrapper:
        "0.74.0/bio/snpeff/download"

# Rule is not submitted as a job to the cluster.
localrules: snpeff_db

# Get the path to the snpeff database. Usually, this is the path to a database from the available
# ones of snpEff. If however a custom db path is given, we return this instead.
def get_snpeff_db_path():
    if config["params"]["snpeff"]["custom-db-dir"]:
        return config["params"]["snpeff"]["custom-db-dir"]
    else:
        return get_snpeff_db_download_path() + config["params"]["snpeff"]["name"]

# =================================================================================================
#     SnpEff
# =================================================================================================

rule snpeff:
    input:
        # (vcf, bcf, or vcf.gz)
        calls="filtered/all.vcf.gz",

        # path to reference db downloaded with the snpeff download wrapper above
        db=get_snpeff_db_path()
    output:
        # annotated calls (vcf, bcf, or vcf.gz)
        calls=report(
            "annotated/snpeff.vcf.gz",
            caption="../reports/vcf.rst",
            category="Calls"
        ),

        # summary statistics (in HTML), optional
        stats=report(
            "annotated/snpeff.html",
            category="Calls"
        ),

        # summary statistics in CSV, optional
        csvstats="annotated/snpeff.csv"
    log:
        "logs/snpeff.log"
    group:
        "snpeff"
    params:
        # optional parameters (e.g., max memory 4g)
        # For finding the chromosome names used by snpeff, add `-v` here
        extra=config["params"]["snpeff"]["extra"]
    conda:
        "../envs/snpeff.yaml"
    wrapper:
        "0.74.0/bio/snpeff/annotate"

# =================================================================================================
#     VEP Downloads
# =================================================================================================

# Get the paths where we store the data, so that we do not have to download it all the time.
# If not provided by the user, we use the directory where the reference genome is.
def get_vep_cache_dir():
    if config["params"]["vep"]["cache-dir"]:
        result = os.path.dirname( config["params"]["vep"]["cache-dir"] )
    else:
        result = os.path.join(
            os.path.dirname( config["data"]["reference-genome"] ), "vep-cache"
        )
    return result

# Same for plugins dir.
def get_vep_plugins_dir():
    if config["params"]["vep"]["plugins-dir"]:
        result = os.path.dirname( config["params"]["vep"]["plugins-dir"] )
    else:
        result = os.path.join(
            os.path.dirname( config["data"]["reference-genome"] ), "vep-plugins"
        )
    return result

rule vep_cache:
    output:
        directory(get_vep_cache_dir()),
    params:
        species  = config["params"]["vep"]["species"],
        release  = config["params"]["vep"]["release"],
        build    = config["params"]["vep"]["build"],
        cacheurl = config["params"]["vep"]["cache-url"],
        # fastaurl = config["params"]["vep"]["fasta-url"],
        # fasta-url: "ftp://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/current/fasta"
    log:
        "logs/vep-cache.log",
    conda:
        # We use a conda environment on top of the wrapper, as the wrapper always causes
        # issues with missing python modules and mismatching program versions and stuff...
        # Here, numpy is missing.
        # Well no, we don't even use the wrapper any more, but our own improved version...
        # See https://github.com/snakemake/snakemake-wrappers/issues/365
        # and https://github.com/snakemake/snakemake-wrappers/issues/366
        "../envs/vep.yaml"
    script:
        "../scripts/vep-cache.py"
    # wrapper:
    #     "0.74.0/bio/vep/cache"

rule vep_plugins:
    output:
        directory(get_vep_plugins_dir()),
    params:
        release = config["params"]["vep"]["release"],
    log:
        "logs/vep-plugins.log",
    conda:
        # Use our own env definition here, to ensure that we are working with the same vep
        # versions across the different rules here. This is not the case in the original wrapper...
        "../envs/vep.yaml"
    wrapper:
        "0.74.0/bio/vep/plugins"

# Rules are not submitted as a job to the cluster.
localrules: vep_cache, vep_plugins

# =================================================================================================
#     VEP
# =================================================================================================

rule vep:
    input:
        calls="filtered/all.vcf.gz",
        cache=get_vep_cache_dir(),
        plugins=get_vep_plugins_dir(),
    output:
        calls=report(
            "annotated/vep.vcf.gz",
            caption="../report/vcf.rst",
            category="Calls",
        ),
        stats=report(
            # The html file has to have a specific file name, so that MultiQC can find it,
            # see https://multiqc.info/docs/#vep
            # At the moment, this is however not yet working, because VEP was only added to
            # MultiQC recently, see https://github.com/ewels/MultiQC/issues/1438
            # Will need to update to MultiQC v1.11 at some point.
            "annotated/vep_summary.html",
            caption="../report/stats.rst",
            category="Calls",
        ),
    params:
        # Pass a list of plugins to use,
        # see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["params"]["vep"]["plugins"],
        extra=config["params"]["vep"]["extra"],
    log:
        "logs/vep-annotate.log",
    threads: 4
    conda:
        # Use our own env definition here, to ensure that we are working with the same vep
        # versions across the different rules here. This is not the case in the original wrapper...
        "../envs/vep.yaml"
    # script:
    #     "../scripts/vep.py"
    # The original snakemake wrapper fails for the Arabidopsis thaliana genome because of an
    # unnecessary/wrong error check, see https://github.com/snakemake/snakemake-wrappers/issues/365
    # It wants there to be one file in the database, when the vep cache wrapper also can produce
    # multiple files/directories in the cache du to also downloading the fasta reference files.
    wrapper:
        "0.74.0/bio/vep/annotate"
