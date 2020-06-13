# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd
import os

# Ensure min Snakemake version
snakemake.utils.min_version("5.7")

# =================================================================================================
#     Basic Configuration
# =================================================================================================

# Add a description of the workflow to the final report
report: "../reports/workflow.rst"

# We check if snakemake was called with what we call the run directory ("rundir"),
# e.g., `snakemake --config rundir="my-run"`, which is the target directory to write all files to.
# If not, we simply write files to the current directory. We explicitly use the run dir instead
# of the working directory offered by snakemake, as the latter makes ALL paths relative to that dir,
# which would mean that we have to re-specify all input file paths as well, and re-load all conda
# modules, etc...
#
# Furthermore, we check if the rundir contains a "config.yaml" configuration file, and load this
# instead of the config.yaml in the main snakemake directory.
# This is useful to have runs that have different settings, but generally re-use the main setup.
if "rundir" in config.keys():
    cfg=os.path.join(config["rundir"], "config.yaml")
    if os.path.isfile(cfg):
        configfile: cfg
    else:
        configfile: "config.yaml"
else:
    config["rundir"] = ""
    configfile: "config.yaml"
if config["rundir"] and not config["rundir"].endswith("/"):
    config["rundir"] += "/"
snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")

# Some settings in the config file need to be converted from empty string to empty list, it seems,
# so that rules that use the files specified in these settings are working properly.
# Maybe there is some better way, but for now, this is working.
if "known-variants" not in config["data"]["reference"] or not config["data"]["reference"]["known-variants"]:
    config["data"]["reference"]["known-variants"]=[]
if "restrict-regions" not in config["settings"] or not config["settings"]["restrict-regions"]:
    config["settings"]["restrict-regions"]=[]

# =================================================================================================
#     Read Samples File
# =================================================================================================

# Read samples and units table
samples = pd.read_table(config["data"]["samples"], dtype=str).set_index(["sample", "unit"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index
snakemake.utils.validate(samples, schema="../schemas/samples.schema.yaml")

# Transform for ease of use
sample_names=list(set(samples.index.get_level_values("sample")))
unit_names=list(set(samples.index.get_level_values("unit")))

# Wildcard constraints: only allow sample names from the spreadsheet to be used
wildcard_constraints:
    sample="|".join(sample_names),
    unit="|".join(unit_names)
    # vartype="snvs|indels" TODO?!

# =================================================================================================
#     Pipeline User Output
# =================================================================================================

# Some helpful messages
logger.info("===========================================================================")
logger.info("    GRENEPIPE")
logger.info("")
logger.info("    Snakefile:          " + (workflow.snakefile))
logger.info("    Base directory:     " + (workflow.basedir))
logger.info("    Working directory:  " + os.getcwd())
logger.info("    Run directory:      " + (config["rundir"][:-1] if config["rundir"] else os.getcwd()))
logger.info("    Config files:       " + (", ".join(workflow.configfiles)))
unitcnt=len(samples.index.get_level_values("unit"))
if unitcnt == len(sample_names):
    logger.info("    Samples:            " + str(len(sample_names)))
else:
    logger.info("    Samples:            " + str(len(sample_names)) + ", with " + str(unitcnt) + " total units")
logger.info("===========================================================================")
logger.info("")

# =================================================================================================
#     Common File Access Functions
# =================================================================================================

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}
