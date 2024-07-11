#!/usr/bin/env python

# =================================================================================================
#     Merge HAFpipe files for all samples and chromosomes
# =================================================================================================

# HAFpipe has the unfortunately oversight in its implementation that its main output files,
# the afSite files with allele frequencies, do not contain any information on their original
# sample or chromosome - except for the file name.
# Hence, our rules expect a certain fixed file naming scheme for these files, and this script
# here expects that same scheme.
# In particular, we expect "{sample}.merged.bam.{chrom}.afSite" for the samples,
# and at least allow to set the base path, which for our case is "hafpipe/allele-frequencies".
# It is hence not terribily portable, but might still be useful for others as inspiration.
# There does not seem to be any other way, given that different use cases might want to handle
# this differently anyway.
# Furthermore, note that we are mostly using shell commands to do the hard work, as this is just
# so much faster than python; we use the snakemake shell command for this - you might need
# to adopt this to some other way of calling if you are working outside of snakemake.

import os
import sys
import math
import resource
import subprocess
from snakemake.shell import shell

# Log everything, and append, to allow us easily to call shell() multiple times
# without having to worry about overwriting the log file.
# We cannot use the snakemake log helper here, as we will change directories, which means
# that the path to the log file will not work. We hence use an absolute path instead.
# log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
if not snakemake.log:
    log = ""
else:
    log = " >> {0} 2>&1".format( os.path.abspath(str(snakemake.log)) )

# -------------------------------------------------------------------------
#     File access
# -------------------------------------------------------------------------

# We actually change to the directory of the files, to keep it simple,
# and to potentially avoid issues with too long paths when doing our hundreds of concats.
os.chdir( snakemake.params.get("base_path", "hafpipe/frequencies") )

# Helper function to get the file path for a HAFpipe afSite file.
def get_afsite_path( sample, chrom ):
    if sample not in snakemake.params.samples:
        raise Exception( "Invalid sample " + sample )
    if chrom not in snakemake.params.chroms:
        raise Exception( "Invalid chrom " + chrom )
    return sample + ".merged.bam." + chrom + ".afSite"

# -------------------------------------------------------------------------
#     Set up and checks
# -------------------------------------------------------------------------

# Get the limit of concurrent files. This is equivalent to `ulimit -n` on Linux.
# If the rule provides an override for this, use that value instead.
soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
if int(snakemake.params.get("concurrent_files", 0)) > 0:
    max_files = int(snakemake.params.get("concurrent_files"))
else:
    # Bit lower than the actual limit, to have room for our output files etc.
    max_files = int(soft) - 10
    if max_files < 500:
        # We probably should do some proper check here if the resource limit retrievel workd...
        # But for now, this is good enough, quick fix to a number that likely works on any system.
        max_files = 500

# Some stupid assertions, just in case, so that we can rely on this downstream.
# Cannot really happen (I think), as the rules would never be executed, but better safe than sorry.
if len(snakemake.params.samples) < 1:
    raise Exception( "No samples provided for HAFpipe merging" )
if len(snakemake.params.chroms) < 1:
    raise Exception( "No chromosomes provided for HAFpipe merging" )

# We are combining all files into larger tables first, and then combine these tables again.
# But we do not do that recursively, so we might hit a limit if there are more of these combined
# tables than we can process in their merging step...
# For a typical system with a `ulimit -n` of 1024, this will happen at ~1mio samples.
# batch_cnt = int(( len(snakemake.params.samples) - 1 ) / max_files) + 1
batch_cnt = int(math.ceil( float(len(snakemake.params.samples)) / float(max_files) ))
if batch_cnt >= max_files:
    raise Exception(
        "Congratulations! You are using incredibly many files (or a very small batch size)! "
        "Unfortunately, at this time, the HAFpipe afSite file merge step in grenepipe "
        "is not implemented for this many files. "
        "Please submit an issue to https://github.com/moiexpositoalonsolab/grenepipe/issues, "
        "so that we know this happend, and can work on a fix for this. Sorry for any inconvenience."
    )

# First, we make sure that all files for one chromosome have the same number of rows,
# otherwise, merging them won't work - we don't do a smart position-based matching here,
# as it does not seem necessary with HAFpipe. If that is wrong, and HAFpipe sometimes validly
# can output files of different length for the same chromosome, we will have to rewrite this.
# This implicitly assumes that the actual positions in these files also match - we do not check
# this as of now, as it seems to be the case, and would be unnecessary effort.
# See https://github.com/petrov-lab/HAFpipe-line/blob/master/calculate_HAFs.sh for evidence.
for chrom in snakemake.params.chroms:
    lc = 0
    for sample in snakemake.params.samples:
        fn=get_afsite_path(sample, chrom)

        # We use the full speed of unix here, instead of the slow python way...
        lines=int(subprocess.check_output("/usr/bin/wc -l " + fn, shell=True).split()[0])
        if lc == 0:
            lc = lines
        if lc != lines:
            raise Exception(
                "Cannot merge HAFpipe afSite files for chromosome " + chrom +
                ", as the per-sample files have different number of positions."
            )

# Some log output
shell(
    "echo -e \"In `pwd` \\n\" {log} ; "
    "echo -e \"Started `date` \\n\" {log} ; "
    "echo \"Merging HAFpipe afSite files\" {log} ; "
    "echo \"Samples: {snakemake.params.samples}\" {log} ; "
    "echo \"Chromosomes: {snakemake.params.chroms}\" {log} ; "
    "echo -e \"Batches: {batch_cnt} \\n\" {log} ; "
)
# print(snakemake.params.samples)
# print(snakemake.params.chroms)

# -------------------------------------------------------------------------
#     Merge chromosomes
# -------------------------------------------------------------------------

# We run the merging first per chromosome, and then simply cat everything together.
# To achieve this, we ditch the header lines, as we will need new ones containing the sample names
# anyway, and we also only use the position column once, and ditch it from all files.
# Potential refactor for the future: We _could_ run this step independently per chromosome,
# in separate executions of a per-chrom rule, and then cat the result in a second rule.
# But given that this is a highly IO intense process, we probably would just spam the disk anyway,
# probably even on clusters. So might not be worth it.
# Also, we always go via the intermediate merged files first, even if potentially the number
# of files is small enough to be merged at once, but then we'd have to write that code twice...
# And for datasets that are that small, the doubling of file copy time probably does not matter
# much anyway, so we safe that effort for now.

chrom_files = ""
for chrom in snakemake.params.chroms:
    # Make a deep copy of the samples list, so that we can delete the ones that we have processed
    # on this chromosome from it... Python and its shallow references...
    samples=list(snakemake.params.samples)

    # Get the first column of the first file, without the header line.
    # This will be the positions column that we use for all, assuming that they all have
    # the same positions (see above).
    sample_0 = get_afsite_path( samples[0], chrom )
    pos_file = "../all-" + chrom + ".pos"
    shell( "cut -d\",\" -f1 {sample_0} | tail -n +2 > {pos_file}" )
    lines=int(subprocess.check_output("/usr/bin/wc -l " + pos_file, shell=True).split()[0])

    # We also need to build a file of the same length that over and over repeats the chromosome,
    # so that we can paste it in front of the chromosome table.
    chr_file = "../all-" + chrom + ".chrom"
    shell(
        "touch {chr_file} ; "
        "for l in `seq 1 {lines}` ; do echo {chrom} >> {chr_file}; done"
    )

    # Some nice user output
    shell(
        "echo \"Processing chromosome {chrom} with {lines} positions\" {log} ; "
    )

    # We process in batches of our max_files size, so that we never have (much) more files open
    # at the same time that than. The output of course as well, but that should work.
    batch_num = 0
    batch_files = ""
    while len(samples) > 0:
        # Move the first max_files many samples into a new list
        batch_smps, samples = samples[:max_files], samples[max_files:]
        batch_smps_cnt = len(batch_smps)

        # We stream each file, using pipe redirections, while at the same time filtering:
        # For each sample, get the frequency column without the head line.
        # We here build a string with all these instructions, and then exectute it.
        sample_paste_files=""
        for sample in batch_smps:
            sn = get_afsite_path( sample, chrom )
            sample_paste_files += " <( cut -d\",\" -f2 " + sn + " | tail -n +2 )"

        # Get just the af column without header of each of the files in the batch.
        batch_file = "../all-" + chrom + "-" + str(batch_num) + ".af"
        shell(
            "echo \" - Batch {batch_num} with {batch_smps_cnt} samples\" {log} ; "
            "paste -d, {sample_paste_files} > {batch_file} ; "
        )
        batch_num += 1
        batch_files += " " + batch_file

    # Now we have the batch files, and need to merge them, including the pos column.
    chrom_file = "../all-" + chrom + ".af"
    shell(
        "echo \" - Merging {batch_num} batch(es)\" {log} ; "
        "paste -d, {chr_file} {pos_file} {batch_files} > {chrom_file} ; "
        "rm {chr_file} {pos_file} {batch_files} ; "
    )
    chrom_files += " " + chrom_file

# -------------------------------------------------------------------------
#     Merge all
# -------------------------------------------------------------------------

# Now we have merged files for each chromosome, all with pos columns, but no headers.
# We now build the header, using our sample names, and then merge the final table.
# We name the columsn after the sample names, and append `.af` to it to mark it as a frequency.

# If we want to compress the table, just add a pipe here. Luckily, appending gzip files to each
# other is a valid thing to do with the format, so we can just build the file chrom by chrom.
if snakemake.params.get("compress", False):
    gzip=" | gzip "
else:
    gzip=""

# Make a header line with all sample names.
header="chrom,pos"
for sample in snakemake.params.samples:
    header += "," + sample + ".af"

# Write the table by pasting all chromsome tables.
all_file = "../all.csv" + (".gz" if gzip else "")
shell(
    "echo \"Merging final table\" {log} ; "
    "echo {header} {gzip} > {all_file} ; "
    "cat {chrom_files} {gzip} >> {all_file} ; "
    "rm {chrom_files} ; "
    "echo -e \"\\nFinished `date`\" {log} ; "
)
