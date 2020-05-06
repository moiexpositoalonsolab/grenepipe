include: "rules/common.smk"

# =================================================================================================
#     All Target Rule
# =================================================================================================

rule all:
    input:
        "genotyped/all.vcf.gz",
        "qc/multiqc.html"

# =================================================================================================
#     Rule Modules
# =================================================================================================

# include: "rules/prep.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
# include: "rules/filtering.smk"
# include: "rules/stats.smk"
include: "rules/qc.smk"
# include: "rules/annotation.smk"
