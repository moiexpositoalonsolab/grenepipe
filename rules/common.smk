# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd
import os
import socket, platform

# Ensure min Snakemake version
snakemake.utils.min_version("5.7")
basedir = workflow.basedir

# =================================================================================================
#     Basic Configuration
# =================================================================================================

# Add a description of the workflow to the final report
report: os.path.join(workflow.basedir, "reports/workflow.rst")

# Load the config. If --directory was provided, this is also loaded from there.
# This is useful to have runs that have different settings, but generally re-use the main setup.
configfile: "config.yaml"
snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")

# Some settings in the config file need to be converted from empty string to empty list, it seems,
# so that rules that use the files specified in these settings are working properly.
# Maybe there is some better way, but for now, this is working.
if "known-variants" not in config["data"]["reference"] or not config["data"]["reference"]["known-variants"]:
    config["data"]["reference"]["known-variants"]=[]
if "restrict-regions" not in config["settings"] or not config["settings"]["restrict-regions"]:
    config["settings"]["restrict-regions"]=[]

# We need to clean up the file name for the reference genome.
if config["data"]["reference"]["genome"].endswith(".gz"):
    config["data"]["reference"]["genome"] = os.path.splitext(config["data"]["reference"]["genome"])[0]

# GATK only accepts reference genomes if their file ending is `fa` or `fasta`, at least in the
# version that we are using. This has been updated later to also include `fas`, see here:
# https://github.com/samtools/htsjdk/commit/657b0a6076d84582a19b741fc28cd3c9a12384bf#diff-77d63fdf4a920459b3a44ead94ad979b93115eba0749baa10a694d061b9d6c1f
# It would of course be better to check the file contents instead of the extenion, but okay...
# So here, we check this, in order to provide some better error messages for users,
# instead of having them run into cyrptic log messages "File is not a supported reference file type".
# Add `".fas", ".fas.gz"` later if we upgrade GATK.
fastaexts = ( ".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna" )
if not config["data"]["reference"]["genome"].endswith( fastaexts ):
    raise Exception(
        "Reference genome file path does not end in " + str(fastaexts) + ", which unfortunately " +
        "is needed by GATK. Please rename the file and change the path in the config.yaml"
    )

# =================================================================================================
#     Read Samples File
# =================================================================================================

# Read samples and units table
samples = pd.read_csv(config["data"]["samples"], sep='\t', dtype=str).set_index(["sample", "unit"], drop=False)
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

# Get a nicely formatted hostname
hostname = socket.gethostname()
hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")

# Some helpful messages
logger.info("===========================================================================")
logger.info("    GRENEPIPE")
logger.info("")
logger.info("    Host:               " + hostname)
logger.info("    Snakefile:          " + (workflow.snakefile))
logger.info("    Base directory:     " + (workflow.basedir))
logger.info("    Working directory:  " + os.getcwd())
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
