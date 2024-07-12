# =================================================================================================
#     Common File Access Functions
# =================================================================================================


# Get the fastq files for a sample, either single or paired end, as a dictionary.
def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    smps = config["global"]["samples"]
    fastqs = smps.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}


# Determine whether a sample is single or paired end.
# We use args to be able to call this function from functions that contain more wildcards
# than just sample and unit, such as the bwa aln rules.
def is_single_end(sample, unit, **kargs):
    """Return True if sample-unit is single end."""
    return pd.isnull(config["global"]["samples"].loc[(sample, unit), "fq2"])


# =================================================================================================
#     Trimming
# =================================================================================================

# Switch to the chosen mapper
trimming_tool_good = False
if config["settings"]["trimming-tool"] == "adapterremoval":

    # Use `adapterremoval`
    include: "trimming-adapterremoval.smk"

    trimming_tool_good = True

elif config["settings"]["trimming-tool"] == "cutadapt":

    # Use `cutadapt`
    include: "trimming-cutadapt.smk"

    trimming_tool_good = True

elif config["settings"]["trimming-tool"] == "fastp":

    # Use `fastp`
    include: "trimming-fastp.smk"

    trimming_tool_good = True

elif config["settings"]["trimming-tool"] == "seqprep":

    # Use `seqprep`
    include: "trimming-seqprep.smk"

    trimming_tool_good = True

elif config["settings"]["trimming-tool"] == "skewer":

    # Use `skewer`
    include: "trimming-skewer.smk"

    trimming_tool_good = True

elif config["settings"]["trimming-tool"] == "trimmomatic":

    # Use `trimmomatic`
    include: "trimming-trimmomatic.smk"

    trimming_tool_good = True

elif config["settings"]["trimming-tool"] == "none":

    # Use dummy implementation that just returns the original fastq files again
    include: "trimming-none.smk"

    trimming_tool_good = True


# Another ugly workaround for https://github.com/snakemake/snakefmt/issues/239
if not trimming_tool_good:
    raise Exception("Unknown trimming-tool: " + config["settings"]["trimming-tool"])
