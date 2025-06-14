__author__ = "Johannes Köster, Lucas Czech"
__copyright__ = "Copyright 2023, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

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

if snakemake.params.get("full_cache_url", ""):
    # Parameters from Snakemake
    url       = snakemake.params.full_cache_url
    outdir    = Path(snakemake.output[0]).expanduser()
    log_file  = snakemake.log[0]

    # Ensure directories exist
    outdir.mkdir(parents=True, exist_ok=True)
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)

    # Derive the tarball filename from the URL
    tarball = outdir / url.rstrip('/').split('/')[-1]

    # Download and extract using Snakemake's shell helper
    shell(
        (
            f"curl -fSL {url} -o {tarball} >> {log_file} 2>&1 && "
            f"tar xzf {tarball} -C {outdir} >> {log_file} 2>&1"
        )
    )

else:
    # Get params. By default, we run only cache (--AUTO c), unlike the original wrapper,
    # which also requestd fasta (--AUTO cf), which would then mess up the check that the
    # subdirectory of the cache contains a single directory that is done in the vep annotation wrapper.
    # See https://github.com/snakemake/snakemake-wrappers/issues/365
    automode = snakemake.params.get("automode", "c")
    extra = snakemake.params.get("extra", "")

    # Extra optional cache and fasta url
    cache_url = snakemake.params.get("cache_url", "")
    if cache_url:
        cache_url = '--CACHEURL "{}"'.format(cache_url)
    fastaurl = snakemake.params.get("fastaurl", "")
    if fastaurl:
        fastaurl = '--FASTAURL "{}"'.format(fastaurl)

    try:
        release = int(snakemake.params.release)
        if snakemake.params.get("cache_release", 0) > 0:
            release = int(snakemake.params.get("cache_release"))
    except ValueError:
        raise ValueError("The parameter release is supposed to be an integer.")

    log = snakemake.log_fmt_shell(stdout=True, stderr=True)

    # Compared to the original wrapper, we add the two urls, and also use a newer version
    # of vep install, which uses --CACHE_VERSION instead of --VERSION.
    shell(
        "vep_install --AUTO {automode} "
        "--SPECIES {snakemake.params.species} "
        "--ASSEMBLY {snakemake.params.build} "
        "--CACHE_VERSION {release} "
        "--CACHEDIR {snakemake.output[0]} "
        "--CONVERT "
        "--NO_UPDATE "
        "{cache_url} {fastaurl} "
        "{extra} {log}"
    )


# import tempfile
# from pathlib import Path
# from snakemake.shell import shell

# # Get params. By default, we run only cache (--AUTO c), unlike the original wrapper,
# # which also requestd fasta (--AUTO cf), which would then mess up the check that the
# # subdirectory of the cache contains a single directory that is done in the vep annotation wrapper.
# # See https://github.com/snakemake/snakemake-wrappers/issues/365
# automode = snakemake.params.get("automode", "c")
# extra = snakemake.params.get("extra", "")
# log = snakemake.log_fmt_shell(stdout=True, stderr=True)


# try:
#     release = int(snakemake.params.release)
#     if snakemake.params.get("cache_release", 0) > 0:
#         release = int(snakemake.params.get("cache_release"))
# except ValueError:
#     raise ValueError("The parameter release is supposed to be an integer.")


# with tempfile.TemporaryDirectory() as tmpdir:
#     # We download the cache tarball manually because vep_install does not consider proxy settings (in contrast to curl).
#     # See https://github.com/bcbio/bcbio-nextgen/issues/1080
#     cache_url = snakemake.params.get("cacheurl", "")
#     if not cache_url:
#         cache_url = snakemake.params.get("url", "ftp://ftp.ensembl.org/pub")
#     cache_tarball = (
#         f"{snakemake.params.species}_vep_{release}_{snakemake.params.build}.tar.gz"
#     )
#     if snakemake.params.get("indexed"):
#         vep_dir = "indexed_vep_cache"
#         convert = ""
#     else:
#         vep_dir = "vep" if snakemake.params.get("url") or release >= 97 else "VEP"
#         convert = "--CONVERT "
#     shell(
#         "curl -L {cache_url}/release-{release}/variation/{vep_dir}/{cache_tarball} -o {tmpdir}/{cache_tarball} {log}"
#     )

#     fastaurl = snakemake.params.get("fastaurl", "")
#     if fastaurl:
#         fastaurl = '--FASTAURL "{}"'.format(fastaurl)

#     log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
#     shell(
#         "vep_install --AUTO {automode} "
#         "--SPECIES {snakemake.params.species} "
#         "--ASSEMBLY {snakemake.params.build} "
#         "--CACHE_VERSION {release} "
#         "--CACHEURL {tmpdir} "
#         "--CACHEDIR {snakemake.output} "
#         "{convert}"
#         "--NO_UPDATE "
#         "{fastaurl} "
#         "{extra} {log}"
#     )