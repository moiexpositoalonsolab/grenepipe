#!/usr/bin/env python
"""
Snakemake SLURM submit script.
"""
import warnings  # use warnings.warn() rather than print() to output info in this script
import os
import datetime
import socket

from snakemake.utils import read_job_properties
import slurm_utils

workingdir = os.getcwd()

# cookiecutter arguments
SBATCH_DEFAULTS = """"""
ADVANCED_ARGUMENT_CONVERSION = {"yes": True, "no": False}["no"]

# Try to find a cluster config for the current host. If not found, reset to empty config file.
CLUSTER_CONFIG = os.path.join(workingdir, "profiles/host/" + socket.gethostname().split('.', 1)[0] + ".yaml")
if not os.path.exists(CLUSTER_CONFIG):
    CLUSTER_CONFIG = ""

RESOURCE_MAPPING = {
    "time": ("time", "runtime", "walltime"),
    "mem": ("mem", "mem_mb", "ram", "memory"),
    "mem-per-cpu": ("mem-per-cpu", "mem_per_cpu", "mem_per_thread"),
    "nodes": ("nodes", "nnodes")
}

# parse job
jobscript = slurm_utils.parse_jobscript()
job_properties = read_job_properties(jobscript)

sbatch_options = {}
cluster_config = slurm_utils.load_cluster_config(CLUSTER_CONFIG)

# 1) sbatch default arguments
sbatch_options.update(slurm_utils.parse_sbatch_defaults(SBATCH_DEFAULTS))

# 2) cluster_config defaults
sbatch_options.update(cluster_config["__default__"])

# 3) Convert resources (no unit conversion!) and threads
sbatch_options.update(
    slurm_utils.convert_job_properties(job_properties, RESOURCE_MAPPING)
)

# 4) cluster_config for particular rule
sbatch_options.update(cluster_config.get(job_properties.get("rule"), {}))

# 5) cluster_config options
sbatch_options.update(job_properties.get("cluster", {}))

# 6) Advanced conversion of parameters
if ADVANCED_ARGUMENT_CONVERSION:
    sbatch_options = slurm_utils.advanced_argument_conversion(sbatch_options)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# Additional features that we want and need.
# Inspiration from: https://github.com/bnprks/snakemake-slurm-profile

# Prepare directory and params for slurm log files
extra_params = {}
extra_params["job_name"] = "snakemake." + job_properties["rule"]
extra_params["log_base"] = os.path.join(workingdir, "slurm-logs")
extra_params["log_dir"] = os.path.join(workingdir, "slurm-logs", job_properties["rule"])
os.makedirs(submission_params["log_dir"], exist_ok=True)

# Set out and err slurm log files
sbatch_options["output"] = "{log_dir}/{job_name}.%j.out".format(**extra_params)
sbatch_options["error"] = "{log_dir}/{job_name}.%j.err".format(**extra_params)

# Write out the submission string. We do not want to rewire too much of this script as of now,
# so instead we do something ugly and duplicate the code from slurm_utils.submit_job() to re-create
# the submission string here. Can beautify in the future.
with open( os.path.join(extra_params["log_base"], "slurm-submissions.log"), "a") as slurmlog:
    now = datetime.datetime.now()
    opt = [f"--{k}={v}" for k, v in sbatch_options.items()]
    cmd = ["sbatch"] + opt + [jobscript]
    slurmlog.write(now.strftime("%Y-%m-%d %H:%M:%S") + "\t" + cmd + "\n")

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#7) Format pattern in snakemake style
sbatch_options = slurm_utils.format_values(sbatch_options, job_properties)

# ensure sbatch output dirs exist
for o in ("output", "error"):
    slurm_utils.ensure_dirs_exist(sbatch_options[o]) if o in sbatch_options else None

# submit job and echo id back to Snakemake (must be the only stdout)
print(slurm_utils.submit_job(jobscript, **sbatch_options))
