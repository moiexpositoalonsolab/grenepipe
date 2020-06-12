# =================================================================================================
#     Dependencies
# =================================================================================================

import os
from snakemake.shell import shell

# =================================================================================================
#     Process Data
# =================================================================================================

known = snakemake.input.get("known", "")
if known:
    known = "--dbsnp " + known

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts", "")
bams = snakemake.input.bam
if isinstance(bams, str):
    bams = [bams]
bams = list(map("-I {}".format, bams))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "gatk --java-options '{java_opts}' HaplotypeCaller {extra} "
    "-R {snakemake.input.ref} {bams} "
    "-ERC GVCF "
    "-O {snakemake.output.gvcf} {known} {log}"
)
