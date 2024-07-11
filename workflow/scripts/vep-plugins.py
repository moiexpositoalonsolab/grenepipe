# =================================================================================================
#     README
# =================================================================================================

# Snakemake wrapper for vep cache, adapted from the wrapper script at
# https://snakemake-wrappers.readthedocs.io/en/0.74.0/wrappers/vep/plugins.html in version 0.74.0
#
# For some non-reproducible reason, we encountered an error with the function outdir.mkdir()
# used in the script: It seems that Snakemake already creates the output directory,
# but that function then fails, as it already exists.
# So here, we use it with `exist_ok = True` to avoid that error...

# =================================================================================================
#     VEP Cache
# =================================================================================================

import sys
from pathlib import Path
from urllib.request import urlretrieve
from zipfile import ZipFile
from tempfile import NamedTemporaryFile

if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

outdir = Path(snakemake.output[0])
outdir.mkdir(exist_ok = True)

with NamedTemporaryFile() as tmp:
    urlretrieve(
        "https://github.com/Ensembl/VEP_plugins/archive/release/{release}.zip".format(
            release=snakemake.params.release
        ),
        tmp.name,
    )

    with ZipFile(tmp.name) as f:
        for member in f.infolist():
            memberpath = Path(member.filename)
            if len(memberpath.parts) == 1:
                # skip root dir
                continue
            targetpath = outdir / memberpath.relative_to(memberpath.parts[0])
            if member.is_dir():
                targetpath.mkdir()
            else:
                with open(targetpath, "wb") as out:
                    out.write(f.read(member.filename))
