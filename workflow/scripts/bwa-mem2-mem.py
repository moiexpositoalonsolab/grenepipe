# =================================================================================================
#     README
# =================================================================================================

# Snakemake wrapper for BWA MEM2, adapted from the wrapper script at
# https://snakemake-wrappers.readthedocs.io/en/v1.7.0/wrappers/bwa-mem2/mem.html,
# but edited, because the original does not incorporate a temporary directory...
# The Snakemake wrapper for bwa mem does, but this one does not :-(
# So we have to take care of this ourselves here, using inspiration from bwa mem,
# see https://snakemake-wrappers.readthedocs.io/en/0.80.0/wrappers/bwa/mem.html

# Also, the original wrapper mixes up idx and index... so buggy...

# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

__author__ = "Christopher Schröder, Johannes Köster, Julian de Ruiter, Lucas Czech"
__copyright__ = (
    "Copyright 2020, Christopher Schröder, Johannes Köster and Julian de Ruiter, 2022,  Lucas Czech"
)
__email__ = (
    "christopher.schroeder@tu-dortmund.de koester@jimmy.harvard.edu, julianderuiter@gmail.com"
)
__license__ = "MIT"

from os import path
import re
import tempfile
from snakemake.shell import shell

# Extract arguments.
extra = " -t {} ".format(snakemake.threads) + snakemake.params.get("extra", "")
sort = snakemake.params.get("sort", "none")
sort_order = snakemake.params.get("sort_order", "coordinate")
sort_extra = snakemake.params.get("sort_extra", "")

# Added from bwa mem wrapper
if re.search(r"-T\b", sort_extra) or re.search(r"--TMP_DIR\b", sort_extra):
    sys.exit(
        "You have specified temp dir (`-T` or `--TMP_DIR`) in params.sort_extra; this is automatically set from params.tmp_dir."
    )

# Also from the bwa mem wrapepr: make sure that the tmp dir is used.
tmp_dir = snakemake.params.get("tmp_dir")
if tmp_dir:
    tempfile.tempdir = tmp_dir

# Prepare index
index = snakemake.input.get("idx", "")
if isinstance(index, str):
    index = path.splitext(snakemake.input.idx)[0]
else:
    index = path.splitext(snakemake.input.idx[0])[0]

# Check inputs/arguments.
if not isinstance(snakemake.input.reads, str) and len(snakemake.input.reads) not in {
    1,
    2,
}:
    raise ValueError("input must have 1 (single-end) or 2 (paired-end) elements")

if sort_order not in {"coordinate", "queryname"}:
    raise ValueError("Unexpected value for sort_order ({})".format(sort_order))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# =================================================================================================
#     Prepare sorting
# =================================================================================================

# Determine which pipe command to use for converting to bam or sorting.
if sort == "none":

    # Simply convert to bam using samtools view.
    pipe_cmd = "samtools view -Sbh -o {snakemake.output[0]} -"

elif sort == "samtools":

    # Sort alignments using samtools sort.
    pipe_cmd = "samtools sort {sort_extra} -T {tmp} -o {snakemake.output[0]} -"

    # Add name flag if needed.
    if sort_order == "queryname":
        sort_extra += " -n"

    # Below are some old fragments to get the tmp dir to properly work.
    # We initially replaced the old prefix temp dir by a proper temp dir.
    # With the original implementation, we would always use the same tmp dir, which causes
    # samtools sort to fail if files from a previous run are still in there...
    # However, we did that as part of sort_extra, introducing a nested {} brace substitution,
    # which is not done by the snakemake shell command, hence using a literal `{tmp}` as the
    # temp dir... so instead, we directly insert this into the command above now.

    # Use prefix for temp.
    # prefix = path.splitext(snakemake.output[0])[0]
    # sort_extra += " -T " + prefix + ".tmp"
    # sort_extra += " -T {tmp}"

elif sort == "picard":

    # Sort alignments using picard SortSam.
    pipe_cmd = (
        "picard SortSam {sort_extra} INPUT=/dev/stdin"
        " OUTPUT={snakemake.output[0]} SORT_ORDER={sort_order}"
    )

else:
    raise ValueError("Unexpected value for params.sort ({})".format(sort))

# =================================================================================================
#     Main Call
# =================================================================================================

# Prepare the shell command
shell_cmd = "(bwa-mem2 mem {extra} {index} {snakemake.input.reads} | " + pipe_cmd + ") {log}"

# When using samtools, we create a temp dir.
if sort == "samtools":
    with tempfile.TemporaryDirectory() as tmp:
        shell(shell_cmd)
else:
    shell(shell_cmd)
