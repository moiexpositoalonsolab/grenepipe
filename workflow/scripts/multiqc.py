# Copied from https://github.com/snakemake/snakemake-wrappers/blob/master/bio/multiqc/wrapper.py,
# respectively (for the exact commit version): https://github.com/snakemake/snakemake-wrappers/blob/d0c7e6894ee33e8bda90e03b157654ff51838876/bio/multiqc/wrapper.py

from os import path
from snakemake.shell import shell

input_dirs = set(path.dirname(fp) for fp in snakemake.input)
output_dir = path.dirname(snakemake.output[0])
output_name = path.basename(snakemake.output[0])
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "multiqc"
    " {snakemake.params}"
    " --force"
    " -o {output_dir}"
    " -n {output_name}"
    " {input_dirs}"
    " {log}"
)
