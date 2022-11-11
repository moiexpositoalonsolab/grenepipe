#!/usr/bin/env python

# =================================================================================================
#     Wrapper script to run HAFpipe in Snakemake
# =================================================================================================

# https://github.com/petrov-lab/HAFpipe-line

import os
import sys
from snakemake.shell import shell

# Log everything, and append, to allow us easily to call shell() multiple times
# without having to worry about overwriting the log file.
# This catches the std out and err, but we also need to provide a log file,
# as HAFpipe prints is normal log to there, see shell command below.
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

# -------------------------------------------------------------------------
#     Check tool availability
# -------------------------------------------------------------------------

# We need a bit of hackery to call HAFpipe, as neither it nor harp (its dependency)
# are available via conda. So we "install" them locally in our grenepipe directory,
# under `grenepipe/packages`. There is an `setup-hafpipe` script there that takes care of this,
# but we need to find this first.
# For all these steps, we get the directory where this python script file is in,
# and then work our way up from there through our grenepipe directory.
script_path = os.path.dirname(os.path.realpath(__file__))
packages_path = os.path.join( script_path, "../packages" )

harp_path   = os.path.join( packages_path, "harp/bin/" )
harp_bin    = os.path.join( packages_path, "harp/bin/harp" )
hafpipe_bin = os.path.join( packages_path, "hafpipe/HAFpipe_wrapper.sh" )

if not os.path.exists( harp_bin ) or not os.path.exists( hafpipe_bin ):
    print("Cannot find harp or HAFpipe in the grenepipe directory. Downloading them now.")
    shell(
        "{packages_path}/setup-hafpipe.sh {log}"
    )

if not os.path.exists( harp_bin ) or not os.path.exists( hafpipe_bin ):
    print(
        "Still cannot find harp or HAFpipe in the grenepipe directory. "
        "Please see grenepipe/packages/README.md to download them manually."
    )
    sys.exit(1)

# -------------------------------------------------------------------------
#     Set up the arguments
# -------------------------------------------------------------------------

# We try to get all hafpipe params that our snakemake rule use, and build the command from that.
args = ""

# `tasks` is always provided via the params, so that we know what this script is supposed to do.
if not snakemake.params.get("tasks", ""):
    raise Exception("Need to provide params: tasks")
args += " --tasks " + snakemake.params.get("tasks")

# `snptable` can be either an output (Task 1), or an input (other Tasks),
# so we check either here, to make sure that we find the right file name.
if snakemake.input.get("snptable", ""):
    args += " --snptable {snakemake.input.snptable:q}"
elif snakemake.output.get("snptable", ""):
    args += " --snptable {snakemake.output.snptable:q}"
else:
    raise Exception("Need to provide input:/output: snptable")

# All other arguments that we use in our rules are taken if provided, for inputs and params.
# We simply list all command line arguments that HAFpipe offers, and check if we find them
# in the places where we expect them from our rules. Outputs are named implicitly in HAFpipe,
# using the inputs as base names, so we do not have to provide any of those, but simply
# need to match them in our rule. Might want to change later, and rename the files to our liking.

# Add potential inputs
hafpipe_inputs = [ "vcf", "bamfile", "refseq" ]
for arg in hafpipe_inputs:
    if snakemake.input.get(arg, ""):
        args += " --" + arg + " {snakemake.input." + arg + ":q}"

# Add potential params, without the `:q`, as those are not files.
hafpipe_params = [ "chrom", "impmethod", "outdir" ]
for arg in hafpipe_params:
    if snakemake.params.get(arg, ""):
        args += " --" + arg + " {snakemake.params." + arg + "}"

# Check that impmethod is actually valid for use in HAFpipe.
# Should not happen with our rules, as they should only call this script for the two valid methods.
impmethod = snakemake.params.get("impmethod", "")
if impmethod != "" and impmethod not in ["simpute", "npute"]:
    raise Exception( "Cannot use impmethod '" + impmethod + "' with HAFpipe directly" )

# Arguments of HAFpipe that exist, but that we do not need to use above:
# --logfile (directly provided below in the shell command)
# --scriptdir (we already know where the scripts are)
# --keephets --subsetlist (extras for Task 1, via config file)
# --nsites (extra for Task 2, via config file)
# --encoding, --generations, --recombrate, --quantile, --winsize (extras for Tasks 3&4, via config file)

# -------------------------------------------------------------------------
#     Run HAFpipe
# -------------------------------------------------------------------------

# Now run the command using the above args, plus other snakemake dependend args.
# We also export the path to harp first, so that it can found by HAFpipe,
# and we finally also forward any extra params that the user might have provided.
shell(
    "export PATH=$PATH:{harp_path} ; "
    "{hafpipe_bin} " + args + \
    "    --logfile {snakemake.log:q} "
    "    {snakemake.params.extra} {log}"
)

# There is a bug in HAFpipe for Task 2 with simpute, where the output file is not named
# as expected, see https://github.com/petrov-lab/HAFpipe-line/issues/4, so here we catch this
# and manually rename to the expected file for that particular case, so that our rule finds it.
# We fixed that in our fork of HAF-pipe, but keep this check around here for completeness.
if snakemake.params.get("tasks") == "2" and impmethod == "simpute":
    if os.path.exists( snakemake.input.snptable + ".imputed" ):
        shell(
            "mv {snakemake.input.snptable:q}.imputed {snakemake.input.snptable:q}.simpute"
        )

# The `numeric_SNPtable.R` script is being run from within Task 1 and produces the `.numeric` file
# of the SNP table. It can require an _insane_ amount of memory for larger founder VCF files,
# but might fail silently when going out of memory.
# In that case, HAF-pipe just continues as if nothing happend, producing a valid, yet empty,
# compressed `.numeric.bgz` file, which will then lead to errors down the line.
# Here, we check this, and issue a proper error for the user, so that they don't have to debug
# this issue, and can just increase their memory.
# We fixed that in our fork of HAF-pipe, but keep this check around here for completeness.
if snakemake.params.get("tasks") == "1":
    # We check that both files exist and that the bgzipped one is not so small that it's likely
    # just a zipped empty file. Reading bgzip in python is tricky, so we don't do that as of now...
    numeric = snakemake.output.get("snptable") + ".numeric"
    numbgz  = numeric + ".bgz"
    if (
        not os.path.exists( numeric ) or
        not os.path.exists( numbgz ) or
        os.path.getsize( numbgz ) < 100
    ):
        raise Exception(
            "The HAF-pipe Task 1 step to convert the SNP table file to a numeric format failed. "
            "It is likely that this is caused by an out-of-memory (OOM) error, as the HAF-pipe "
            "R script `numeric_SNPtable.R` reads in the whole SNP table at once. This can be quite "
            "big if your founder VCF has many samples. Please check all log files for errors, "
            "and try to increase the amount of memory for our grenepipe rule `hafpipe_snp_table`."
        )
