# =================================================================================================
#     Dependencies
# =================================================================================================

import pandas as pd

# =================================================================================================
#     General Snakemake Pipeline Settings
# =================================================================================================

# Ensure min Snakemake version
snakemake.utils.min_version("4.3.1")

# TODO
# report: "../report/workflow.rst"

# Get input config and sample files, and validate them
configfile: "config.yaml"
# TODO
# snakemake.utils.validate(config, schema="../schemas/config.schema.yaml")

# Set optional settings of config to usable python values.
# Otherwise, we have to do this all over the place...
if "known-variants" not in config["data"]["reference"] or not config["data"]["reference"]["known-variants"]:
    config["data"]["reference"]["known-variants"]=[]
if "restrict-regions" not in config["data"]["reference"] or not config["data"]["reference"]["restrict-regions"]:
    config["data"]["reference"]["restrict-regions"]=[]

# =================================================================================================
#     Read Samples File
# =================================================================================================

# Read samples and units table
samples = pd.read_table(config["data"]["samples"], dtype=str).set_index(["sample", "unit"], drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])  # enforce str in index
snakemake.utils.validate(samples, schema="../schemas/samples.schema.yaml")

# Transform for ease of use
sample_indices=(samples.index.get_level_values("sample"))
unit_indices=(samples.index.get_level_values("unit"))

# Wildcard constraints: only allow sample names from the spreadsheet to be used
wildcard_constraints:
    sample="|".join(sample_indices),
    unit="|".join(unit_indices)
    # vartype="snvs|indels" TODO?!

# =================================================================================================
#     Common File Access Functions
# =================================================================================================

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}

# =================================================================================================
#     Common Param Functions
# =================================================================================================

def get_regions_param(regions=config["data"]["reference"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["settings"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default

def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, default="--intervals {}".format(wildcards.contig)) +
            config["params"]["gatk"]["HaplotypeCaller"])
