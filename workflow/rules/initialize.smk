# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd
import os, sys, pwd, re
import socket, platform
import subprocess
from datetime import datetime
import logging
from pathlib import Path

from snakemake_interface_executor_plugins.settings import ExecMode

# Ensure min Snakemake version
snakemake.utils.min_version("8.15.2")
basedir = workflow.basedir

# We are currently setting up our own extra log file, so that the below banner is shown.
# Snakemake currently only activates logging to the `.snakemake/log` files _after_ having
# processed all snakefiles, which is not really how logging should work...
# See https://github.com/snakemake/snakemake/issues/2974 for the issue.
# We need to distinguish between the main instance, and the instances of each rule job.
if logger.mode == ExecMode.DEFAULT:
    extra_logdir = "snakemake"
else:
    extra_logdir = "snakemake-jobs"
os.makedirs(os.path.join("logs", extra_logdir), exist_ok=True)
extra_logfile = os.path.abspath(
    os.path.join(
        "logs",
        extra_logdir,
        datetime.now().isoformat().replace(":", "") + ".log",
    )
)
logger.logger.addHandler(logging.FileHandler(extra_logfile))


# =================================================================================================
#     Basic Configuration
# =================================================================================================

# We want to report the grenepipe version for the user, for reproducibility.
# The following line is automatically replaced by the deploy scripts. Do not change manually.
grenepipe_version = "0.13.4"  # GRENEPIPE_VERSION #


# Load the config. If --directory was provided, this is also loaded from there.
# This is useful to have runs that have different settings, but generally re-use the main setup.
configfile: "config.yaml"


# After changing to our new scheme, we can verify the scheme to fit our expextation.
snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")


# Add a description of the workflow to the final report
report: os.path.join(workflow.basedir, "reports/workflow.rst")


# Include the functions neeed to initialize the pipeline for analysing a set of fastq samples,
# or, if the mappings table is given, starting from there.
# This reads the samples/mappings table, and provides validation and user output functions for it.
include: "initialize-reference.smk"


if "mappings-table" in config["data"] and config["data"]["mappings-table"]:

    include: "initialize-bam.smk"

else:

    include: "initialize-fastq.smk"


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
logger.info("    Samples:            " + get_sample_units_print())
logger.info("")
logger.info("=====================================================================================")
logger.info("")

# No need to have these output vars available in the rest of the snakefiles
del indent
del pltfrm, hostname, username
del conda_ver, conda_env
del cmdline, cfgfiles
