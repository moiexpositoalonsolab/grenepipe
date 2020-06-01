include: "rules/common.smk"

# =================================================================================================
#     Default "All" Target Rule
# =================================================================================================

# The rule that is executed by default. We include the result files of different intermediate steps
# here as well, for example, the genotyped vcf file from the "calling.smk" step, so that a nice
# arrow shows up in the DAG that reminds us that this is an important intermediate file.
rule all:
    input:
        # Basic steps
        config["rundir"] + "genotyped/all.vcf.gz",
        config["rundir"] + "filtered/all.vcf.gz",
        config["rundir"] + "annotated/all.vcf.gz",

        # Quality control
        config["rundir"] + "qc/multiqc.html",
        config["rundir"] + "plots/depths.svg",
        config["rundir"] + "plots/allele-freqs.svg"

# =================================================================================================
#     Rule Modules
# =================================================================================================

# The preparation rule is special (for now), and needs to be run beforehand.
# See the file itself for instructions on how to run it.
# include: "rules/prep.smk"

include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/annotation.smk"
include: "rules/qc.smk"
include: "rules/stats.smk"
