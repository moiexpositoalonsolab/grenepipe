include: "rules/common.smk"

# =================================================================================================
#     All Target Rule
# =================================================================================================

# The rule that is executed by default. We include the result files of different steps here as well,
# for example, the genotyped vcf file from the "calling.smk" step, so that it shows up in the DAG.
rule all:
    input:
        config["rundir"] + "genotyped/all.vcf.gz",
        config["rundir"] + "filtered/all.vcf.gz",
        config["rundir"] + "qc/multiqc.html"

# =================================================================================================
#     Rule Modules
# =================================================================================================

# The preparation rule is special (for now), and needs to be run beforehand.
# See the file itself for instructions on how to run it.
# include: "rules/prep.smk"

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
# include: "rules/stats.smk"
include: "rules/qc.smk"
# include: "rules/annotation.smk"
