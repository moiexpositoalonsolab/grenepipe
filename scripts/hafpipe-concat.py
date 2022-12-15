#!/usr/bin/env python

# =================================================================================================
#     Concatenate HAFpipe files for one sample
# =================================================================================================

# This script concatenates the per-chromosome afSite files of one sample into one file,
# containing columns for the chromosome, position, and frequency.
# It hence is yet another mitigation of a shortcoming of HAFpipe, which only produces
# single files per chromosome, with no information on which chromosome it is in the file itself...

import os
import sys
import math
import resource
import subprocess
from snakemake.shell import shell

# -------------------------------------------------------------------------
#     Setup
# -------------------------------------------------------------------------

# Log everything, and append, to allow us easily to call shell() multiple times
# without having to worry about overwriting the log file.
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

# User output.
shell(
    "echo \"Started `date`\" {log} ; "
    "echo \"Sample: {snakemake.params.sample}\" {log} ; "
    "echo \"Chromosomes: {snakemake.params.chroms}\" {log} ; "
)

# Helper function to get the file path for a HAFpipe afSite file.
def get_afsite_path( sample, chrom ):
    if chrom not in snakemake.params.chroms:
        raise Exception( "Invalid chrom " + chrom )
    return "hafpipe/frequencies/" + sample + ".bam." + chrom + ".afSite"

# -------------------------------------------------------------------------
#     Concat chromosomes
# -------------------------------------------------------------------------

# If we want to compress the table, just add a pipe here. Luckily, appending gzip files to each
# other is a valid thing to do with the format, so we can just build the file chrom by chrom.
if snakemake.params.get("compress", False):
    gzip=" | gzip "
else:
    gzip=""

# Make a header line with the sample name, and write it to the output.
sample_file = snakemake.output.table
header="chrom,pos," + snakemake.params.sample + ".af"
shell(
    "echo {header} {gzip} > {sample_file} ; "
)

# Go through all chromosomes in the specified order, and concat them to the output file,
# while prepending the chromosome as a first column to the file, and leaving out its header.
for chrom in snakemake.params.chroms:
    sn = get_afsite_path( snakemake.params.sample, chrom )
    shell(
        "cat {sn} | tail -n +2 | sed 's/^/{chrom},/' {gzip} >> {sample_file} ; "
    )

# User output.
shell(
    "echo \"Finished `date`\" {log} ; "
)
