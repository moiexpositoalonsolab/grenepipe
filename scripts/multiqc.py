# =================================================================================================
#     Dependencies
# =================================================================================================

from os import path
from snakemake.shell import shell

# =================================================================================================
#     Process Data
# =================================================================================================

input_dirs = set(path.dirname(fp) for fp in snakemake.input)
output_dir = path.dirname(snakemake.output[0])
output_name = path.basename(snakemake.output[0])
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "export LC_ALL=C.UTF-8; "
    "export LANG=C.UTF-8; "
    "multiqc {snakemake.params} --force -o {output_dir} -n {output_name} {input_dirs} {log}"
)
