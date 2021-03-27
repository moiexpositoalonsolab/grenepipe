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
        "genotyped/all.vcf.gz",
        "filtered/all.vcf.gz",
        "annotated/all.vcf.gz" if config["settings"]["snpeff"] else [],

        # Quality control
        "qc/multiqc.html",

        # Stats. Some deactivated for now, as they run out of memory for our dataset
        # "plots/depths.svg",
        # "plots/allele-freqs.svg",
        # "tables/sample-sizes.tsv"
        "tables/frequencies.tsv" if config["settings"]["frequency-table"] else []

# The main `all` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all

# =================================================================================================
#     All QC, but not SNP calling
# =================================================================================================

# This alternative target rule executes all quality control (QC) steps of read trimming and mapping,
# but does not call SNPs, and does not call snpeff. The result is mainly the MultiQC report (without
# the snpeff part however), as well as the fastqc reports.
rule all_qc:
    input:
        # Quality control
        "qc/multiqc.html"

# The `all_qc` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all_qc

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
include: "rules/damage.smk"
