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
    "if (locale -a | grep \"C.UTF-8\") >/dev/null 2>&1 ; then      \n"
    "    export LC_ALL=C.UTF-8                                     \n"
    "    export LANG=C.UTF-8                                       \n"
    "elif (locale -a | grep \"en_US.utf8\") >/dev/null 2>&1 ; then \n"
    "    export LC_ALL=en_US.utf8                                  \n"
    "    export LANG=en_US.utf8                                    \n"
    "fi                                                            \n"
    "multiqc {snakemake.params} --force -o {output_dir} -n {output_name} {input_dirs} {log}"
)
