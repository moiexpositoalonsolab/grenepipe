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

# Prepare directory and params for slurm log files
workingdir = os.getcwd()
extra_params = {}
extra_params["log_base"] = os.path.join(workingdir, "slurm-logs")
os.makedirs(extra_params["log_base"], exist_ok=True)

def write_debug_log(msg):
    with open( os.path.join(extra_params["log_base"], "slurm-debug.log"), "a") as debug:
        now = datetime.datetime.now()
        debug.write(now.strftime("%Y-%m-%d %H:%M:%S") + "\t" + str(msg) + "\n")
    # pass

# cookiecutter arguments
SBATCH_DEFAULTS = """"""
ADVANCED_ARGUMENT_CONVERSION = {"yes": True, "no": False}["no"]

# Try to find a cluster config in the dir where the current script is located at;
# this supposedly also works if we use a symlink to the script.
# If no config file found, reset to empty config file.
CLUSTER_CONFIG = os.path.join( os.path.dirname(__file__), "cluster_config.yaml" )
if not os.path.exists(CLUSTER_CONFIG):
    CLUSTER_CONFIG = ""

RESOURCE_MAPPING = {
    "time": ("time", "runtime", "walltime", "time_min"),
    "mem": ("mem", "mem_mb", "ram", "memory"),
    "mem-per-cpu": ("mem-per-cpu", "mem_per_cpu", "mem_per_thread"),
    "nodes": ("nodes", "nnodes")
}

# parse job
jobscript = slurm_utils.parse_jobscript()
job_properties = read_job_properties(jobscript)

sbatch_options = {}
cluster_config = slurm_utils.load_cluster_config(CLUSTER_CONFIG)
# write_debug_log( "c\t" + str(cluster_config))

# 1) sbatch default arguments
sbatch_options.update(slurm_utils.parse_sbatch_defaults(SBATCH_DEFAULTS))
write_debug_log( "1\t" + str(sbatch_options))

# 2) cluster_config defaults
sbatch_options.update(cluster_config["__default__"])
write_debug_log( "2\t" + str(sbatch_options))

# 3) Convert resources (no unit conversion!) and threads
sbatch_options.update(
    slurm_utils.convert_job_properties(job_properties, RESOURCE_MAPPING)
)
write_debug_log( "3\t" + str(sbatch_options))

# 4) cluster_config for particular rule or group
if job_properties["type"] == "single":
    sbatch_options.update(cluster_config.get(job_properties.get("rule"), {}))
elif job_properties["type"] == "group":
    sbatch_options.update(cluster_config.get(job_properties.get("groupid"), {}))
else:
    print("Error: slurm-submit.py doesn't support job type {} yet!".format(job_properties["type"]))
    sys.exit(1)
write_debug_log( "4\t" + str(sbatch_options))

# 5) cluster_config options
sbatch_options.update(job_properties.get("cluster", {}))
write_debug_log( "5\t" + str(sbatch_options))

# 6) Advanced conversion of parameters
if ADVANCED_ARGUMENT_CONVERSION:
    sbatch_options = slurm_utils.advanced_argument_conversion(sbatch_options)
write_debug_log( "6\t" + str(sbatch_options))

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# Additional features that we want and need.
# Inspiration from: https://github.com/bnprks/snakemake-slurm-profile

def file_escape(string):
    return string.replace("/", "_").replace(" ", "_")

# Prepare job name for log script
if job_properties["type"] == "single":
    extra_params["job_name"] = "snakejob." + job_properties["rule"]
    extra_params["log_dir"] = os.path.join(workingdir, "slurm-logs", job_properties["rule"])
elif job_properties["type"] == "group":
    extra_params["job_name"] = "snakejob." + job_properties["groupid"]
    extra_params["log_dir"] = os.path.join(workingdir, "slurm-logs", job_properties["groupid"])
else:
    print("Error: slurm-submit.py doesn't support job type {} yet!".format(job_properties["type"]))
    sys.exit(1)
if "wildcards" in job_properties and len(job_properties["wildcards"]) > 0:
    extra_params["job_name"] += "." + ".".join([key + "=" + file_escape(value) for key,value in job_properties["wildcards"].items()])
os.makedirs(extra_params["log_dir"], exist_ok=True)

# Set job name and out and err slurm log files
sbatch_options["job-name"] = extra_params["job_name"]
sbatch_options["output"] = "{log_dir}/{job_name}.%j.out".format(**extra_params)
sbatch_options["error"] = "{log_dir}/{job_name}.%j.err".format(**extra_params)
write_debug_log( "S\t" + str(sbatch_options) + "\n")

# Write out the submission string. We do not want to rewire too much of this script as of now,
# so instead we do something ugly and duplicate the code from slurm_utils.submit_job() to re-create
# the submission string here. Can beautify in the future.
with open( os.path.join(extra_params["log_base"], "slurm-submissions.log"), "a") as slurmlog:
    now = datetime.datetime.now()
    opt = ["--"+str(k)+"="+str(v) for k, v in sbatch_options.items()]
    cmd = ["sbatch"] + opt + [jobscript]
    slurmlog.write(now.strftime("%Y-%m-%d %H:%M:%S") + "\t" + extra_params["job_name"] + "\t" + ' '.join(cmd) + "\n")

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#7) Format pattern in snakemake style
sbatch_options = slurm_utils.format_values(sbatch_options, job_properties)

# ensure sbatch output dirs exist
for o in ("output", "error"):
    slurm_utils.ensure_dirs_exist(sbatch_options[o]) if o in sbatch_options else None

# submit job and echo id back to Snakemake (must be the only stdout)
print(slurm_utils.submit_job(jobscript, **sbatch_options))
