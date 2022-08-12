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

# These rules need to come after the `all` rule above, as Snakemake takes the first rule it finds
# as the implicit target when no other target is specified, and we want that to be the `all` rule.

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
include: "rules/frequency.smk"
