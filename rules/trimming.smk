# =================================================================================================
#     Trimming
# =================================================================================================

# Switch to the chosen mapper
if config["settings"]["trimming-tool"] == "cutadapt":

    # Use `cutadapt`
    include: "trimming-cutadapt.smk"

elif config["settings"]["trimming-tool"] == "trimmomatic":

    # Use `trimmomatic`
    include: "trimming-trimmomatic.smk"

elif config["settings"]["trimming-tool"] == "skewer":

    # Use `skewer`
    include: "trimming-skewer.smk"

else:
    raise Exception("Unknown trimming-tool: " + config["settings"]["trimming-tool"])
