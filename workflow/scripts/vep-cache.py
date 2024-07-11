# =================================================================================================
#     README
# =================================================================================================

# Snakemake wrapper for vep cache, adapted from the wrapper script at
# https://github.com/snakemake/snakemake-wrappers/blob/master/bio/vep/cache/wrapper.py in version 0.74.0,
# https://github.com/snakemake/snakemake-wrappers/blob/6b71f64fba7ee2c6cad31315d9ccb1ed26c4605c/bio/vep/cache/wrapper.py
#
# We here fix some issues with the original script before the wrapper is fixed.
# In particular, we make the requirement that fasta is downloaded optional,
# and add the capability to set download URLs for the cache and fasta files.

# =================================================================================================
#     VEP Cache
# =================================================================================================

from pathlib import Path
from snakemake.shell import shell

# Get params. By default, we run only cache (--AUTO c), unlike the original wrapper,
# which also requestd fasta (--AUTO cf), which would then mess up the check that the
# subdirectory of the cache contains a single directory that is done in the vep annotation wrapper.
# See https://github.com/snakemake/snakemake-wrappers/issues/365
automode = snakemake.params.get("automode", "c")
extra = snakemake.params.get("extra", "")

# Extra optional cache and fasta url
cacheurl = snakemake.params.get("cacheurl", "")
if cacheurl:
    cacheurl = "--CACHEURL \"{}\"".format(cacheurl)
fastaurl = snakemake.params.get("fastaurl", "")
if fastaurl:
    fastaurl = "--FASTAURL \"{}\"".format(fastaurl)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Compared to the original wrapper, we add the two urls, and also use a newer version
# of vep install, which uses --CACHE_VERSION instead of --VERSION.
shell(
    "vep_install --AUTO {automode} "
    "--SPECIES {snakemake.params.species} "
    "--ASSEMBLY {snakemake.params.build} "
    "--CACHE_VERSION {snakemake.params.release} "
    "--CACHEDIR {snakemake.output} "
    "--CONVERT "
    "--NO_UPDATE "
    "{cacheurl} {fastaurl} "
    "{extra} {log}"
)
