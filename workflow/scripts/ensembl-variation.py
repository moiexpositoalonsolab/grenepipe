# =================================================================================================
#     README
# =================================================================================================

# Snakemake wrapper for ensembl-variation, adapted from the wrapper script at
# https://github.com/snakemake/snakemake-wrappers/blob/master/bio/reference/ensembl-variation/wrapper.py
# in version "v3.13.6/bio/reference/ensembl-variation"
#
# This script serves as a fix for https://github.com/snakemake/snakemake-wrappers/issues/3070
# Instead of trying to resolve the ensembl path from its parts, we just download the file directly.
# Way more stable.

# =================================================================================================
#     ensembl-variation
# =================================================================================================

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import tempfile
import subprocess
import sys
import os
from snakemake.shell import shell
from snakemake.exceptions import WorkflowError

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Get the download url and name
url = snakemake.params.get("url")
name = os.path.basename(url)

try:
    shell('echo "Trying download from {url}" {log} ; ')
    gather = "curl -O {url}".format(url=url)
    workdir = os.getcwd()
    out = os.path.abspath(snakemake.output[0])
    with tempfile.TemporaryDirectory() as tmpdir:
        if snakemake.input.get("fai"):
            fai = os.path.abspath(snakemake.input.fai)
            shell(
                "(cd {tmpdir}; {gather} && "
                "bcftools concat -Oz --naive {name} > concat.vcf.gz && "
                "bcftools reheader --fai {fai} concat.vcf.gz "
                "> {out}) {log}"
            )
        else:
            shell(
                "(cd {tmpdir}; {gather} && "
                "bcftools concat -Oz --naive {name} "
                "> {out}) {log}"
            )
except subprocess.CalledProcessError as e:
    if snakemake.log:
        sys.stderr = open(snakemake.log[0], "a")
    print(
        "Unable to download variation data from Ensembl.",
        file=sys.stderr,
    )
    exit(1)
