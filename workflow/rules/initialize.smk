# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd
import os, sys, pwd, re
import socket, platform
import subprocess
from datetime import datetime

# Ensure min Snakemake version
snakemake.utils.min_version("5.7")
basedir = workflow.basedir

# =================================================================================================
#     Basic Configuration
# =================================================================================================

# We want to report the grenepipe version for the user, for reproducibility.
# The following line is automatically replaced by the deploy scripts. Do not change manually.
grenepipe_version = "0.12.2"  # GRENEPIPE_VERSION#


# Add a description of the workflow to the final report
report: os.path.join(workflow.basedir, "reports/workflow.rst")


# Load the config. If --directory was provided, this is also loaded from there.
# This is useful to have runs that have different settings, but generally re-use the main setup.
configfile: "config.yaml"


# Legacy fixes to avoid disrupting user workflows that already have a config file from a previous
# version of grenepipe. We change the config keys here to our new scheme.
# We renamed "samples" to "samples-table":
if "samples" in config["data"] and "samples-table" not in config["data"]:
    config["data"]["samples-table"] = config["data"]["samples"]
    del config["data"]["samples"]

# After changing to our new scheme, we can verify the scheme to fit our expextation.
snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")

# Some settings in the config file need to be converted from empty string to empty list, it seems,
# so that rules that use the files specified in these settings are working properly.
# Maybe there is some better way, but for now, this is working.
# We need to do this after the scheme validation, as we are introducing changes here
# that are not according to the scheme.
if "known-variants" not in config["data"] or not config["data"]["known-variants"]:
    config["data"]["known-variants"] = []
if "restrict-regions" not in config["settings"] or not config["settings"]["restrict-regions"]:
    config["settings"]["restrict-regions"] = []

# We need to clean up the file name for the reference genome.
# The prep rule decompress_genome provides the unzipped genome as needed.
if config["data"]["reference-genome"].endswith(".gz"):
    config["data"]["reference-genome"] = os.path.splitext(config["data"]["reference-genome"])[0]

# GATK only accepts reference genomes if their file ending is `fa` or `fasta`, at least in the
# version that we are using. This has been updated later to also include `fas`, see here:
# https://github.com/samtools/htsjdk/commit/657b0a6076d84582a19b741fc28cd3c9a12384bf#diff-77d63fdf4a920459b3a44ead94ad979b93115eba0749baa10a694d061b9d6c1f
# It would of course be better to check the file contents instead of the extension, but okay...
# So here, we check this, in order to provide some better error messages for users,
# instead of having them run into cryptic log messages "File is not a supported reference file type".
# Add `".fas", ".fas.gz"` later if we upgrade GATK.
fastaexts = (".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fas", ".fas.gz")
if not config["data"]["reference-genome"].endswith(fastaexts):
    raise Exception(
        "Reference genome file path does not end in "
        + str(fastaexts)
        + ", which unfortunately is needed by GATK to be able to find the file. "
        + "Please rename the file and change the path in the config.yaml"
    )

# =================================================================================================
#     Read Samples Table
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
    config["data"]["samples-table"], sep="\t", dtype=str
).set_index(["sample", "unit"], drop=False)
config["global"]["samples"].index = config["global"]["samples"].index.set_levels(
    [i.astype(str) for i in config["global"]["samples"].index.levels]
)
snakemake.utils.validate(config["global"]["samples"], schema="../schemas/samples.schema.yaml")


# Helper function to get a list of all units of a given sample name.
def get_sample_units(sample):
    res = list()
    for unit in config["global"]["samples"].loc[sample].unit:
        if unit not in res:
            res.append(unit)
    return res


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
config["global"]["unit-names"] = list(set(config["global"]["samples"].index.get_level_values("unit")))


# Wildcard constraints: only allow sample names from the spreadsheet to be used
wildcard_constraints:
    sample="|".join(config["global"]["sample-names"]),
    unit="|".join(config["global"]["unit-names"]),


# =================================================================================================
#     Pipeline User Output
# =================================================================================================

# The final output is tabular, we might need to indent subsequent lines correctly.
indent = 24

# Get a nicely formatted username and hostname
username = pwd.getpwuid(os.getuid())[0]
hostname = socket.gethostname()
hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")

# Get some info on the platform and OS
pltfrm = platform.platform() + "\n" + (" " * indent) + platform.version()
try:
    # Not available in all versions, so we need to catch this
    ld = platform.linux_distribution()
    if len(ld):
        pltfrm += "\n" + (" " * indent) + ld
    del ld
except:
    pass
try:
    # Mac OS version comes back as a nested tuple?!
    # Need to merge the tuples...
    def merge_tuple(x, bases=(tuple, list)):
        for e in x:
            if type(e) in bases:
                for e in merge_tuple(e, bases):
                    yield e
            else:
                yield e

    mv = " ".join(merge_tuple(platform.mac_ver()))
    if not mv.isspace():
        pltfrm += "\n" + (" " * indent) + mv
    del mv, merge_tuple
except:
    pass

# Get the git commit hash of grenepipe, if available.
try:
    process = subprocess.Popen(
        ["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = process.communicate()
    out = out.decode("ascii")
    grenepipe_git_hash = out.strip()
    if grenepipe_git_hash:
        grenepipe_version += "-" + grenepipe_git_hash
    del process, out, err, grenepipe_git_hash
except:
    pass

# Get the conda version, if available.
try:
    process = subprocess.Popen(["conda", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode("ascii")
    conda_ver = out[out.startswith("conda") and len("conda") :].strip()
    del process, out, err
    if not conda_ver:
        conda_ver = "n/a"
except:
    conda_ver = "n/a"

# Same for mamba. This somehow can also give a differing conda version.
# Who knows what that means. I'm sick of conda. Just reporting the version here,
# and have someone else deal with it.
try:
    process = subprocess.Popen(["mamba", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode("ascii")
    mamba_ver = re.findall("mamba *(.*) *", out)[0]
    conda_ver_mamba = re.findall("conda *(.*) *", out)[0]
    del process, out, err
    if not mamba_ver:
        mamba_ver = "n/a"
        conda_ver_mamba = ""
except:
    mamba_ver = "n/a"
    conda_ver_mamba = ""
if conda_ver_mamba and conda_ver_mamba != conda_ver:
    conda_ver += " (conda), " + conda_ver_mamba + " (mamba)"

# Get the conda env name, if available.
# See https://stackoverflow.com/a/42660674/4184258
conda_env = os.environ["CONDA_DEFAULT_ENV"] + " (" + os.environ["CONDA_PREFIX"] + ")"
if conda_env == " ()":
    conda_env = "n/a"

# Get nicely wrapped command line
cmdline = sys.argv[0]
for i in range(1, len(sys.argv)):
    if sys.argv[i].startswith("--"):
        cmdline += "\n" + (" " * indent) + sys.argv[i]
    else:
        cmdline += " " + sys.argv[i]

# Get abs paths of all config files
cfgfiles = []
for cfg in workflow.configfiles:
    cfgfiles.append(os.path.abspath(cfg))
cfgfiles = "\n                        ".join(cfgfiles)

# Get a nice output of the number of samples and units
unitcnt = len(config["global"]["samples"].index.get_level_values("unit"))
if unitcnt == len(config["global"]["sample-names"]):
    smpcnt = str(len(config["global"]["sample-names"]))
else:
    smpcnt = str(len(config["global"]["sample-names"])) + ", with " + str(unitcnt) + " total units"

# Main grenepipe header, helping with debugging etc for user issues
logger.info("=====================================================================================")
logger.info(r"       _____         _______ __   __   _______ ______  ___   ______   _______ ")
logger.info(r"      /  ___\ ____  /  ____//  \ /  / /  ____/|   _  \ \  \ |   _  \ /  ____/ ")
logger.info(r"     /  /____|  _ \|  |___  |   \|  ||  |___  |  |_]  ||  | |  |_]  |  |___   ")
logger.info(r"    |  /|__  | |_) |   ___| |       ||   ___| |   ___/ |  | |   ___/|   ___|  ")
logger.info(r"    \  \__|  |  _ <|  |____ |  |\   ||  |____ |  |     |  | |  |    |  |____  ")
logger.info(r"     \______/|_| \_\_______\/__| \__|\_______\|__|     \___\|__|    \_______\ ")
logger.info("")
logger.info("    Date:               " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
logger.info("    Platform:           " + pltfrm)
logger.info("    Host:               " + hostname)
logger.info("    User:               " + username)
logger.info("    Conda:              " + str(conda_ver))
logger.info("    Mamba:              " + str(mamba_ver))
logger.info("    Python:             " + str(sys.version.split(" ")[0]))
logger.info("    Snakemake:          " + str(snakemake.__version__))
logger.info("    Grenepipe:          " + str(grenepipe_version))
logger.info("    Conda env:          " + str(conda_env))
logger.info("    Command:            " + cmdline)
logger.info("")
logger.info("    Base directory:     " + workflow.basedir)
logger.info("    Working directory:  " + os.getcwd())
logger.info("    Config file(s):     " + cfgfiles)
logger.info("    Samples:            " + smpcnt)
logger.info("")
logger.info("=====================================================================================")
logger.info("")

# No need to have these output vars available in the rest of the snakefiles
del indent
del pltfrm, hostname, username
del conda_ver, conda_env
del cmdline, cfgfiles
del unitcnt, smpcnt

# =================================================================================================
#     Sample and File Validity Checks
# =================================================================================================

# Check that the number of samples fits with the expectation, if provided by the user.
if config["data"].get("samples-count", 0) > 0:
    if config["data"]["samples-count"] != len(config["global"]["sample-names"]):
        raise Exception(
            "Inconsistent number of samples found in the sample table ("
            + str(len(config["global"]["sample-names"]))
            + ") and the samples count consistency check provided in the config file ("
            + str(config["data"]["samples-count"])
            + ")."
        )

# Check the three file paths provided in the config file directly.
if not os.path.isabs(config["data"]["samples-table"]):
    logger.warning(
        "Path to the samples table as provided in the config file is not an absolute path. "
        "We recommend using absolute paths for all files.\n"
    )
if not os.path.isabs(config["data"]["reference-genome"]):
    logger.warning(
        "Path to the reference genome as provided in the config file is not an absolute path. "
        "We recommend using absolute paths for all files.\n"
    )
if config["data"]["known-variants"] and not os.path.isabs(config["data"]["known-variants"]):
    logger.warning(
        "Path to samples table as provided in the config file is not an absolute path. "
        "We recommend using absolute paths for all files.\n"
    )


# Helper to check that a string contains no invalid chars for file names.
# This is just the file, not its path! Slashes are considered invalid by this function.
def valid_filename(fn):
    # Only accept alnum, underscore, dash, and dot.
    return fn.replace("_", "").replace("-", "").replace(".", "").isalnum() and fn.isascii()
    # Bit more loose: allow all except Windows forbidden chars.
    # return not( True in [c in fn for c in "<>:\"/\\|?*"])


# Helper to check if a file path contains weird characters.
# We just want to warn about this for input fastq files, but still try to continue.
def valid_filepath(fn):
    # Only accept alnum, underscore, dash, dot, and slashes.
    clean = fn.replace("_", "").replace("-", "").replace(".", "").replace("/", "").replace("\\", "")
    return clean.isalnum() and clean.isascii()


# List that contains tuples for all samples with their units.
# In other words, a list of tuples of the sample and unit column of the sample table,
# in the same order.
config["global"]["sample-units"] = list()
problematic_filenames = 0
relative_filenames = 0
for index, row in config["global"]["samples"].iterrows():
    if (row["sample"], row["unit"]) in config["global"]["sample-units"]:
        raise Exception(
            "Multiple rows with identical sample name and unit found in samples table: "
            + str(row["sample"])
            + " "
            + str(row["unit"])
        )
    config["global"]["sample-units"].append((row["sample"], row["unit"]))

    # Do a check that the sample and unit names are valid file names.
    # They are used for file names, and would cause weird errors if they contain bad chars.
    if not valid_filename(row["sample"]) or not valid_filename(row["unit"]):
        raise Exception(
            "Invalid sample name or unit name found in samples table that contains characters "
            + "which cannot be used as sample/unit names for naming output files: "
            + str(row["sample"])
            + " "
            + str(row["unit"])
            + "; for maximum robustness, we only allow alpha-numerical, dots, dashes, and underscores. "
            + "Use for example the script `tools/copy-samples.py --mode link [...] --clean` "
            + "to create a new samples table and symlinks to the existing fastq files to solve this."
        )

    # Do a check of the fastq file names.
    if not os.path.isfile(row["fq1"]) or (not pd.isnull(row["fq2"]) and not os.path.isfile(row["fq2"])):
        raise Exception(
            "Input fastq files listed in the input files table "
            + config["data"]["samples-table"]
            + " not found: "
            + str(row["fq1"])
            + "; "
            + str(row["fq2"])
        )
    if not valid_filepath(row["fq1"]) or (not pd.isnull(row["fq2"]) and not valid_filepath(row["fq2"])):
        problematic_filenames += 1
    if not os.path.isabs(row["fq1"]) or (not pd.isnull(row["fq2"]) and not os.path.isabs(row["fq2"])):
        relative_filenames += 1

# Warning about input names and files.
if problematic_filenames > 0:
    logger.warning(
        str(problematic_filenames)
        + " of the "
        + str(len(config["global"]["sample-names"]))
        + " input samples listed in the input files table "
        + config["data"]["samples-table"]
        + " contain problematic characters. We generally advise to only use alpha-numeric "
        "characters, dots, dashes, and underscores. "
        "Use for example the script `tools/copy-samples.py --mode link [...] --clean` "
        "to create a new samples table and symlinks to the existing fastq files to solve this. "
        "We will try to continue running with these files, but it might lead to errors.\n"
    )
if relative_filenames > 0:
    logger.warning(
        str(relative_filenames)
        + " of the "
        + str(len(config["global"]["sample-names"]))
        + " input samples listed in the input files table "
        + config["data"]["samples-table"]
        + " use relative file paths. We generally advise to only use absolute paths. "
        "Use for example the script `tools/generate-table.py` "
        "to create a samples table with absolute paths. "
        "We will try to continue running with these files, but it might lead to errors.\n"
    )

del problematic_filenames, relative_filenames


# Check if a given string can be converted to a number, https://stackoverflow.com/q/354038/4184258
def is_number(s):
    try:
        float(s)  # for int, long, float
    except ValueError:
        return False
    return True


# We check if any sample names are only numbers, and warn about this. Bioinformatics is messy...
numeric_sample_names = 0
for sn in config["global"]["sample-names"]:
    if is_number(sn):
        numeric_sample_names += 1

if numeric_sample_names > 0:
    logger.warning(
        str(numeric_sample_names)
        + " of the "
        + str(len(config["global"]["sample-names"]))
        + " input sample names listed in the input files table "
        + config["data"]["samples-table"]
        + " are just numbers, or can be converted to numbers. We generally advise to avoid this, "
        "as it might confuse downstream processing such as Excel and R. The same applies for names "
        'that can be converted to dates ("Dec21"), but we do not check this here. '
        "We will now continue running with these files, as we can work with this here, "
        "but recommend to change the sample names.\n"
    )

del numeric_sample_names
