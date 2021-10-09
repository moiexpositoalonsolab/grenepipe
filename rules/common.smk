# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd
import os, sys, pwd
import socket, platform
import subprocess
from datetime import datetime

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
# The prep rule decompress_genome provides the unzipped genome as needed.
if config["data"]["reference"]["genome"].endswith(".gz"):
    config["data"]["reference"]["genome"] = os.path.splitext(config["data"]["reference"]["genome"])[0]

# GATK only accepts reference genomes if their file ending is `fa` or `fasta`, at least in the
# version that we are using. This has been updated later to also include `fas`, see here:
# https://github.com/samtools/htsjdk/commit/657b0a6076d84582a19b741fc28cd3c9a12384bf#diff-77d63fdf4a920459b3a44ead94ad979b93115eba0749baa10a694d061b9d6c1f
# It would of course be better to check the file contents instead of the extension, but okay...
# So here, we check this, in order to provide some better error messages for users,
# instead of having them run into cryptic log messages "File is not a supported reference file type".
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

# We add the samples information to the config, in order to not spam our global scope
# (unfortunately, while python is super good with namespaces, it is super bad with scopes,
# and in particular in snakemake, every file is included so that all variables defined in a file
# are global and accessible in all subsequent files as well...).
# We use a new top level key that is not used in the config file for this. Assert this.
if "global" in config:
    raise Exception("Config key 'global' already defined. Someone messed with our setup.")
else:
    config["global"] = {}

# Read samples and units table, and enforce to use strings in the index
config["global"]["samples"] = pd.read_csv(
    config["data"]["samples"], sep='\t', dtype=str).set_index(["sample", "unit"], drop=False
)
config["global"]["samples"].index = config["global"]["samples"].index.set_levels(
    [i.astype(str) for i in config["global"]["samples"].index.levels]
)
snakemake.utils.validate( config["global"]["samples"], schema="../schemas/samples.schema.yaml" )

# Get a list of all samples names, in the same order as the input sample table.
# Samples with multiple units appear only once, at the first position in the table.
# We cannot use a simple approach here, as this messes up the sample
# order, which we do not want... (good that we noticed that bug though!)
# So instead, we iterate, and add sample names incrementally.
config["global"]["sample-names"] = list()
for index, row in config["global"]["samples"].iterrows():
    s = row["sample"]
    if s not in config["global"]["sample-names"]:
        config["global"]["sample-names"].append(s)

# Unordered list of all unit names that appear in all the samples.
config["global"]["unit-names"] = list(set(
    config["global"]["samples"].index.get_level_values("unit")
))

# List that contains tuples for all samples with their units.
# In other words, a list of tuples of the sample and unit column of the sample table,
# in the same order.
config["global"]["sample-units"] = list()
for index, row in config["global"]["samples"].iterrows():
    if ( row["sample"], row["unit"] ) in config["global"]["sample-units"]:
        raise Exception(
            "Identical sample name and unit found in samples table: " +
            str(row["sample"]) + " " + str(row["unit"])
        )
    config["global"]["sample-units"].append(( row["sample"], row["unit"] ))

# Helper function to get a list of all units of a given sample name.
def get_sample_units( sample ):
    res = list()
    for unit in config["global"]["samples"].loc[sample].unit:
        if unit not in res:
            res.append(unit)
    return res

# Wildcard constraints: only allow sample names from the spreadsheet to be used
wildcard_constraints:
    sample="|".join(config["global"]["sample-names"]),
    unit="|".join(config["global"]["unit-names"])
    # vartype="snvs|indels" TODO?!

# =================================================================================================
#     Pipeline User Output
# =================================================================================================

# Get a nicely formatted username and hostname
username = pwd.getpwuid( os.getuid() )[ 0 ]
hostname = socket.gethostname()
hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")

# Get the conda version, if available.
try:
    process = subprocess.Popen(['conda', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode('ascii')
    conda_ver = out[out.startswith("conda") and len("conda"):].strip()
    del process, out, err
    if not conda_ver:
        conda_ver = "n/a"
except:
    conda_ver = "n/a"

# Get nicely wrapped command line
cmdline = sys.argv[0]
for i in range( 1, len(sys.argv)):
    if sys.argv[i].startswith("--"):
        cmdline += "\n                        " + sys.argv[i]
    else:
        cmdline += " " + sys.argv[i]

# Get abs paths of all config files
cfgfiles = []
for cfg in workflow.configfiles:
    cfgfiles.append( os.path.abspath(cfg) )
cfgfiles = ", ".join(cfgfiles)

# Get a nice output of the number of samples and units
unitcnt=len(config["global"]["samples"].index.get_level_values("unit"))
if unitcnt == len(config["global"]["sample-names"]):
    smpcnt = str(len(config["global"]["sample-names"]))
else:
    smpcnt = str(len(config["global"]["sample-names"])) + ", with " + str(unitcnt) + " total units"

# Some helpful messages
logger.info("=====================================================================================")
logger.info("    GRENEPIPE")
logger.info("")
logger.info("    Date:               " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
logger.info("    User:               " + username)
logger.info("    Host:               " + hostname)
logger.info("    Conda:              " + str(conda_ver))
logger.info("    Python:             " + str(sys.version.split(' ')[0]))
logger.info("    Snakemake:          " + str(snakemake.__version__))
logger.info("    Command:            " + cmdline)
logger.info("")
logger.info("    Base directory:     " + workflow.basedir)
logger.info("    Working directory:  " + os.getcwd())
logger.info("    Config file(s):     " + cfgfiles)
logger.info("    Samples:            " + smpcnt)
logger.info("=====================================================================================")
logger.info("")

# No need to have these output vars available in the rest of the snakefiles
del username
del hostname
del conda_ver
del cmdline
del cfgfiles
del unitcnt
del smpcnt

# =================================================================================================
#     Common File Access Functions
# =================================================================================================

# Get the fastq files for a sample, either single or paired end, as a dictionary.
def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = config["global"]["samples"].loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}

# Determine whether a sample is single or paired end.
# We use args to be able to call this function from functions that contain more wildcards
# than just sample and unit, such as the bwa aln rules.
def is_single_end( sample, unit, **kargs ):
    """Return True if sample-unit is single end."""
    return pd.isnull(config["global"]["samples"].loc[(sample, unit), "fq2"])
