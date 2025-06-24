# =================================================================================================
#     Dependencies
# =================================================================================================

from datetime import datetime
from pathlib import Path
import inspect
import logging
import os, sys, pwd, re
import pandas as pd
import socket, platform
import subprocess
import yaml

from snakemake_interface_executor_plugins.settings import ExecMode

# Ensure min Snakemake version
snakemake.utils.min_version("8.15.2")
basedir = workflow.basedir

# =================================================================================================
#     Fix Snakemake Logging
# =================================================================================================

# We are currently setting up our own extra log file, so that the below banner is shown.
# Snakemake currently only activates logging to the `.snakemake/log` files _after_ having
# processed all snakefiles, which is not really how logging should work...
# See https://github.com/snakemake/snakemake/issues/2974 for the issue.
# Furthermore, it does not respeced info messages any more in Snakemake 9.3.3,
# see https://github.com/snakemake/snakemake/issues/3558
# Hence, yet again, we need our own tooling to fix other people's mistakes...

# First, set up a custom log file that we can present to our users.
# We need to distinguish between the main instance, and the instances of each rule job.
if logger_manager.mode is None or logger_manager.mode == ExecMode.DEFAULT:
    extra_logdir = "snakemake"
elif logger_manager.mode == ExecMode.SUBPROCESS:
    extra_logdir = "snakemake-subprocess"
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
logger_manager.logger.addHandler(logging.FileHandler(extra_logfile))


# For now, we define our own wrapper around the wrapper of the snakemake logging...
# That allows us to use this as a single point of modification should they finally
# manage to get the logging to work properly.
# Right now as of snakemake 9.3.3, info level is not printed at all...
# See https://github.com/snakemake/snakemake/issues/3558
# So for now, we promote everything to a warning... that sucks, but otherwise,
# our users would not be able to see the grenepipe header etc.
# Furthermore, we print everything to terminal as well, because we have to.
def fix_log_info(message):
    logger.info(message)
    # print(message, file=sys.stdout)
    sys.__stdout__.write(message + "\n")
    sys.__stdout__.flush()


def fix_log_warn(message):
    logger.warning(message)
    # print(message, file=sys.stdout)
    sys.__stdout__.write(message + "\n")
    sys.__stdout__.flush()


# =================================================================================================
#     Basic Configuration
# =================================================================================================

# We want to report the grenepipe version for the user, for reproducibility.
# The following line is automatically replaced by the deploy scripts. Do not change manually.
grenepipe_version = "0.15.0"  # GRENEPIPE_VERSION #


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
#     Resource Configuration
# =================================================================================================


# Obtain the resources.yaml file. First, we check the path specified in the config.yaml.
# If that is empty, we check the working directory. If that also does not contain a resources
# file, we fall back to the default one in the grenepipe directory.
resources_file = config["settings"].get("resources-yaml", "")
if resources_file and not os.path.isfile(resources_file):
    raise Exception("Invalid path to resources.yaml specified in config.yaml: " + resources_file)
if not resources_file:
    if os.path.isfile("resources.yaml"):
        resources_file = "resources.yaml"
    else:
        resources_file = workflow.basedir + "/../config/resources.yaml"
if not resources_file or not os.path.isfile(resources_file):
    raise Exception("Coud not find resources.yaml")

with open(resources_file) as f:
    resources_config = yaml.safe_load(f)


# Helper function to get the number of cpus specified i nthe resource config.
# Unfortunately, we need to set this for every rule, as snakemake processes the threads
# on the first pass already, and so we cannot set it later any more. Any later changes
# would not correctly affect the thread allocations.
# For the mem and time resources, this is different, as those are not first-class resources
# of snakemake, and so we can set them in bulk later (at the end of the main Snakefile).
def get_rule_threads(rule_name):
    default = resources_config["default"]["cpus"]
    return int(resources_config.get(rule_name, {}).get("cpus", default))


# =================================================================================================
#     Pipeline User Output
# =================================================================================================


# The final output is tabular, we might need to indent subsequent lines correctly.
indent = 24

# Get some info on the platform and OS
def info_platform():
    pltfrm = platform.platform() + "\n" + (" " * indent) + platform.version()
    try:
        # Not available in all versions, so we need to catch this
        ld = platform.linux_distribution()
        if len(ld):
            pltfrm += "\n" + (" " * indent) + ld
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
    except:
        pass
    return pltfrm

# Get a nicely formatted hostname
def info_hostname():
    hostname = socket.gethostname()
    hostname = hostname + ("; " + platform.node() if platform.node() != socket.gethostname() else "")
    return hostname

# Get a nicely formatted username
def info_username():
    return pwd.getpwuid(os.getuid())[0]

# Get the conda version, if available.
def info_conda_version():
    conda_ver = "n/a"
    try:
        process = subprocess.Popen(
            ["conda", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = process.communicate()
        out = out.decode("ascii")
        conda_ver = out[out.startswith("conda") and len("conda") :].strip()
        if not conda_ver:
            conda_ver = "n/a"
    except:
        pass
    return str(conda_ver)

# Same for mamba. This somehow can also give a differing conda version.
# Who knows what that means. I'm sick of conda. Just reporting the version here,
# and have someone else deal with it.
def info_mamba_version():
    mamba_ver = ""

    # Get normal mamba first, for the old mamba version output
    try:
        process = subprocess.Popen(
            ["mamba", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = process.communicate()
        out = out.decode("ascii")
        mamba_ver = re.findall("mamba *(.*) *", out)[0]
        conda_ver_mamba = re.findall("conda *(.*) *", out)[0]
        if not mamba_ver:
            mamba_ver = ""
            conda_ver_mamba = ""
    except:
        mamba_ver = ""
        conda_ver_mamba = ""
    if conda_ver_mamba and conda_ver_mamba != info_conda_version():
        mamba_ver += ", with conda " + conda_ver_mamba

    # If that did not work, try the new mamba version output
    if not mamba_ver:
        try:
            process = subprocess.Popen(
                ["mamba", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            out, err = process.communicate()
            mamba_ver = str(out.decode("ascii")).strip()
        except:
            mamba_ver = ""

    # Lastly, also check for micromamba, for full info
    try:
        process = subprocess.Popen(
            ["micromamba", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = process.communicate()
        micromamba_ver = str(out.decode("ascii")).strip()
    except:
        micromamba_ver = ""
    if micromamba_ver:
        if mamba_ver:
            mamba_ver += "\n    Micromamba:         " + micromamba_ver
        else:
            mamba_ver = micromamba_ver + " (micromamba)"

    # Finaly check: if we did not find any mamba, report n/a
    if not mamba_ver:
        mamba_ver = "n/a"
    return str(mamba_ver)

# Get the grenepipe version including git commit hash of grenepipe, if available.
def info_grenepipe_version():
    gpv = grenepipe_version
    try:
        process = subprocess.Popen(
            ["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        out, err = process.communicate()
        out = out.decode("ascii")
        grenepipe_git_hash = out.strip()
        if grenepipe_git_hash:
            gpv += "-" + grenepipe_git_hash
    except:
        pass
    return str(gpv)

# Get the python version currently executing this script.
# If this differs in sub-instances of snakemake, we have an issue.
def info_python_version():
    return str(sys.version.split(" ")[0])

# Get the snakemake version running this script.
def info_snakemake_version():
    return str(snakemake.__version__)

# Get the conda env name, if available.
# See https://stackoverflow.com/a/42660674/4184258
def info_conda_env():
    conda_env = os.environ["CONDA_DEFAULT_ENV"] + " (" + os.environ["CONDA_PREFIX"] + ")"
    if conda_env == " ()":
        conda_env = "n/a"
    return str(conda_env)

# Get nicely wrapped command line
def info_command_line():
    cmdline = sys.argv[0]
    for i in range(1, len(sys.argv)):
        if sys.argv[i].startswith("--"):
            cmdline += "\n" + (" " * indent) + sys.argv[i]
        else:
            cmdline += " " + sys.argv[i]
    return cmdline

# Get abs paths of all config files
def info_config_files():
    cfgfiles = []
    for cfg in workflow.configfiles:
        cfgfiles.append(os.path.abspath(cfg))
    if resources_file:
        cfgfiles.append(os.path.abspath(resources_file))
    cfgfiles = "\n                        ".join(cfgfiles)
    return cfgfiles

# Main grenepipe header, helping with debugging etc for user issues
fix_log_info("=====================================================================================")
fix_log_info(r"       _____         _______ __   __   _______ ______  ___   ______   _______ ")
fix_log_info(r"      /  ___\ ____  /  ____//  \ /  / /  ____/|   _  \ \  \ |   _  \ /  ____/ ")
fix_log_info(r"     /  /____|  _ \|  |___  |   \|  ||  |___  |  |_]  ||  | |  |_]  |  |___   ")
fix_log_info(r"    |  /|__  | |_) |   ___| |       ||   ___| |   ___/ |  | |   ___/|   ___|  ")
fix_log_info(r"    \  \__|  |  _ <|  |____ |  |\   ||  |____ |  |     |  | |  |    |  |____  ")
fix_log_info(r"     \______/|_| \_\_______\/__| \__|\_______\|__|     \___\|__|    \_______\ ")
fix_log_info("")
fix_log_info("    Date:               " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
fix_log_info("    Platform:           " + info_platform())
fix_log_info("    Host:               " + info_hostname())
fix_log_info("    User:               " + info_username())
fix_log_info("    Conda:              " + info_conda_version())
fix_log_info("    Mamba:              " + info_mamba_version())
fix_log_info("    Python:             " + info_python_version())
fix_log_info("    Snakemake:          " + info_snakemake_version())
fix_log_info("    Grenepipe:          " + info_grenepipe_version())
fix_log_info("    Conda env:          " + info_conda_env())
fix_log_info("    Command:            " + info_command_line())
fix_log_info("")
fix_log_info("    Base directory:     " + workflow.basedir)
fix_log_info("    Working directory:  " + os.getcwd())
fix_log_info("    Config file(s):     " + info_config_files())
fix_log_info("    Samples:            " + get_sample_units_print())
fix_log_info("")
fix_log_info("=====================================================================================")
fix_log_info("")


# No need to have these vars available in the rest of the snakefiles
del indent