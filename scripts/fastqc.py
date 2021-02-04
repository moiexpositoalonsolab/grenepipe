# =================================================================================================
#     README
# =================================================================================================

# Snakemake wrapper for fastqc, adapted from the wrapper script at
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastqc.html,
# but edited, because the original just causes too much trouble due to the tmp dir
# (e.g., in cluster environments), and also silences the logging, which is bad for debugging.
# We here do not use the --quite option, so that we can see what went wrong.

# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from os import path
import re
from tempfile import TemporaryDirectory
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

def basename_without_ext(file_path):
    """Returns basename of file path, without the file extension."""

    base = path.basename(file_path)

    # Original line that does not work with fq and the like
    # split_ind = 2 if base.endswith(".fastq.gz") else 1
    # base = ".".join(base.split(".")[:-split_ind])

    # Minimalistic fix to at least get fq extensions to work
    # split_ind = 2 if base.endswith(".fastq.gz") or base.endswith(".fq.gz") else 1
    # base = ".".join(base.split(".")[:-split_ind])

    # Proper solution that follows the actual fastqc code.
    # Proposed for the snakemake wrapper that this script is based on:
    # https://github.com/snakemake/snakemake-wrappers/issues/288
    base = re.sub('\\.gz$', '', base)
    base = re.sub('\\.bz2$', '', base)
    base = re.sub('\\.txt$', '', base)
    base = re.sub('\\.fastq$', '', base)
    base = re.sub('\\.fq$', '', base)
    base = re.sub('\\.sam$', '', base)
    base = re.sub('\\.bam$', '', base)

    return base

# =================================================================================================
#     Main Call
# =================================================================================================

# Run fastqc, since there can be race conditions if multiple jobs
# use the same fastqc dir, we create a temp dir.
# That means we also need to take special care of the log file, as its dir might not be present
# in the temp dir, so also use a dummy here and move later.
with TemporaryDirectory() as tempdir:

    # More verbose output than the wrapper, for debugging
    shell(
        # "echo \"{tempdir:q}\" ; "
        "echo \"Input:       {snakemake.input[0]:q}\" > {tempdir:q}/log.txt ; "
        "echo \"Output zip:  {snakemake.output.zip:q}\" >> {tempdir:q}/log.txt ; "
        "echo \"Output html: {snakemake.output.html:q}\" >> {tempdir:q}/log.txt ; "
        "echo -e \"\n--\n\" >> {tempdir:q}/log.txt ; "
        "fastqc {snakemake.params} -t {snakemake.threads} "
        "--outdir {tempdir:q} {snakemake.input[0]:q} >> {tempdir:q}/log.txt 2>&1 ;"
        # "ls {tempdir:q}"
    )

    # Move outputs into proper position.
    output_base = basename_without_ext(snakemake.input[0])
    html_path = path.join(tempdir, output_base + "_fastqc.html")
    zip_path = path.join(tempdir, output_base + "_fastqc.zip")
    log_path = path.join(tempdir, "log.txt")

    if snakemake.output.html != html_path:
        shell("mv {html_path:q} {snakemake.output.html:q}")

    if snakemake.output.zip != zip_path:
        shell("mv {zip_path:q} {snakemake.output.zip:q}")

    if snakemake.log != log_path:
        shell("mv {log_path:q} {snakemake.log:q}")
