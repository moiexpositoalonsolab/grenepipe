# =================================================================================================
#     Common
# =================================================================================================

# We first need to load the common functionality, which gives us access to the config file,
# generates the sample list, and prepares some other things for us that are needed below.
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
        "annotated/snpeff.vcf.gz" if config["settings"]["snpeff"] else [],
        "annotated/vep.vcf.gz" if config["settings"]["vep"] else [],

        # Quality control
        "qc/multiqc.html",

        # Reference genome statistics
        config["data"]["reference"]["genome"] + ".seqkit",

        # Pileup
        expand(
            "mpileup/{sample}-individual-units.mpileup.gz",
            sample=config["global"]["sample-names"]
        ) if "samples-individual-units" in config["settings"]["pileups"] else [],
        expand(
            "mpileup/{sample}-merged-units.mpileup.gz",
            sample=config["global"]["sample-names"]
        ) if "samples-merged-units" in config["settings"]["pileups"] else [],
        "mpileup/all-individual-units.mpileup.gz" if "all-individual-units" in config["settings"]["pileups"] else [],
        "mpileup/all-merged-units.mpileup.gz" if "all-merged-units" in config["settings"]["pileups"] else [],
        "mpileup/all-merged-samples.mpileup.gz" if "all-merged-samples" in config["settings"]["pileups"] else [],

        # Stats. Some deactivated for now.
        "tables/frequencies.tsv" if config["settings"]["frequency-table"] else []
        # "plots/depths.svg",
        # "plots/allele-freqs.svg",
        # "tables/sample-sizes.tsv"

# The main `all` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all

# =================================================================================================
#     Rule Modules
# =================================================================================================

include: "rules/prep.smk"
include: "rules/trimming.smk"
include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/annotation.smk"
include: "rules/qc.smk"
include: "rules/stats.smk"
include: "rules/damage.smk"
include: "rules/pileup.smk"

# =================================================================================================
#     All QC, but not SNP calling
# =================================================================================================

# This alternative target rule executes all quality control (QC) steps of read trimming and mapping,
# but does not call SNPs, and does not call snpeff. The result is mainly the MultiQC report (without
# the snpeff part however), as well as the fastqc reports.
rule all_qc:
    input:
        # Quality control
        "qc/multiqc.html",

        # Reference genome statistics
        config["data"]["reference"]["genome"] + ".seqkit"

# The `all_qc` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all_qc

# =================================================================================================
#     All bams, but not SNP calling
# =================================================================================================

# This alternative target rule executes all steps up to th mapping, and yields the final bam
# files that would otherwise be used for variant calling in the downstream process.
# That is, depending on the config, these are the sorted, filtered, remove duplicates, or
# recalibrated base qualities bam files.
rule all_bams:
    input:
        get_all_bams()

# The `all_bams` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all_bams

# =================================================================================================
#     All pileups, but not SNP calling
# =================================================================================================

# This alternative target rule executes all steps up to th mapping and mpileup creation.
# This is the same as the above all_bams rule, but additionally also requests the pileups.
rule all_pileups:
    input:
        expand(
            "mpileup/{sample}-individual-units.mpileup.gz",
            sample=config["global"]["sample-names"]
        ) if "samples-individual-units" in config["settings"]["pileups"] else [],
        expand(
            "mpileup/{sample}-merged-units.mpileup.gz",
            sample=config["global"]["sample-names"]
        ) if "samples-merged-units" in config["settings"]["pileups"] else [],
        "mpileup/all-individual-units.mpileup.gz" if "all-individual-units" in config["settings"]["pileups"] else [],
        "mpileup/all-merged-units.mpileup.gz" if "all-merged-units" in config["settings"]["pileups"] else [],
        "mpileup/all-merged-samples.mpileup.gz" if "all-merged-samples" in config["settings"]["pileups"] else [],

# The `all_pileups` rule is local. It does not do anything anyway,
# except requesting the other rules to run.
localrules: all_pileups
