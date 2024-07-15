# =================================================================================================
#     Read Mapping Table
# =================================================================================================

# Similar setup to what we do in `initialize-fastq.smk`, see there for details.
if "global" in config:
    raise Exception("Config key 'global' already defined. Someone messed with our setup.")
else:
    config["global"] = {}

# Read mappings table, and enforce to use strings in the index
config["global"]["samples"] = pd.read_csv(
    config["data"]["mappings-table"], sep="\t", dtype=str
).set_index(["sample"], drop=False)

# Get a list of all samples names, in the same order as the input sample table.
# We cannot use a simple approach here, as this messes up the sample
# order, which we do not want... (good that we noticed that bug though!)
# So instead, we iterate, and add sample names incrementally.
config["global"]["sample-names"] = list()
for index, row in config["global"]["samples"].iterrows():
    s = row["sample"]
    if s not in config["global"]["sample-names"]:
        config["global"]["sample-names"].append(s)


# Helper function for user output to summarize the samples found in the table.
# In the fastq case, we here print the units as well.
# For bam files, it's just simply the number of samples.
def get_sample_units_print():
    return str(len(config["global"]["sample-names"]))


# Wildcard constraints: only allow sample names from the spreadsheet to be used
wildcard_constraints:
    sample="|".join(config["global"]["sample-names"]),


# =================================================================================================
#     Sample and File Validity Checks
# =================================================================================================


# Check that the number of samples fits with the expectation, if provided by the user.
if config["data"].get("samples-count", 0) > 0:
    if config["data"]["samples-count"] != len(config["global"]["sample-names"]):
        raise Exception(
            "Inconsistent number of samples found in the sample table ("
            + str(len(config["global"]["sample-names"]))
            + ") and the samples count consistency check provided in the config file under "
            + "`data: samples-count` ("
            + str(config["data"]["samples-count"])
            + ")."
        )

# We recommend to use absolute paths. Check that for the samples table.
if not os.path.isabs(config["data"]["mappings-table"]):
    logger.warning(
        "Path to the samples table as provided in the config file is not an absolute path. "
        "We recommend using absolute paths for all files.\n"
    )

# Run integrity checks similar to what we do for the fastq files.
uniq_sample_names = list()
problematic_filenames = 0
relative_filenames = 0
for index, row in config["global"]["samples"].iterrows():
    if row["sample"] in uniq_sample_names:
        raise Exception(
            "Multiple rows with identical sample name found in samples table: " + str(row["sample"])
        )
    uniq_sample_names.append(row["sample"])

    # Do a check that the sample names are valid file names.
    # They are used for file names, and would cause weird errors if they contain bad chars.
    if not valid_filename(row["sample"]):
        raise Exception(
            "Invalid sample name name found in samples table that contains characters "
            + "which cannot be used as sample names for naming output files: "
            + str(row["sample"])
            + "; for maximum robustness, we only allow alpha-numerical, dots, dashes, "
            "and underscores. "
            + "Use for example the script `tools/copy-samples.py --mode link [...] --clean` "
            + "to create a new samples table and symlinks to the existing fastq files to solve this."
        )

    # Do a check of the fastq file names.
    if not os.path.isfile(row["bam"]):
        raise Exception(
            "Input bam files listed in the input files mappings table "
            + config["data"]["mappings-table"]
            + " not found: "
            + str(row["bam"])
        )
    if not valid_filepath(row["bam"]):
        problematic_filenames += 1
    if not os.path.isabs(row["bam"]):
        relative_filenames += 1

# Warning about input names and files.
if problematic_filenames > 0:
    logger.warning(
        str(problematic_filenames)
        + " of the "
        + str(len(config["global"]["sample-names"]))
        + " input samples listed in the input files mappings table "
        + config["data"]["mappings-table"]
        + " contain problematic characters. We generally advise to only use alpha-numeric "
        "characters, dots, dashes, and underscores. "
        "We will try to continue running with these files, but it might lead to errors.\n"
    )
if relative_filenames > 0:
    logger.warning(
        str(relative_filenames)
        + " of the "
        + str(len(config["global"]["sample-names"]))
        + " input samples listed in the input files mappings table "
        + config["data"]["mappings-table"]
        + " use relative file paths. We generally advise to only use absolute paths. "
        "We will try to continue running with these files, but it might lead to errors.\n"
    )


# Check if a given string can be converted to a number, https://stackoverflow.com/q/354038/4184258
def is_number(s):
    try:
        float(s)  # for int, long, float
    except ValueError:
        return False
    return True


# We check if any sample names are only numbers, and warn about this. Bioinformatics is messy...
numeric_sample_names = 0
for sn in config["global"]["sample-names"]:
    if is_number(sn):
        numeric_sample_names += 1

if numeric_sample_names > 0:
    logger.warning(
        str(numeric_sample_names)
        + " of the "
        + str(len(config["global"]["sample-names"]))
        + " input sample names listed in the input files mappings table "
        + config["data"]["mappings-table"]
        + " are just numbers, or can be converted to numbers. We generally advise to avoid this, "
        "as it might confuse downstream processing such as Excel and R. The same applies for names "
        'that can be converted to dates ("Dec21"), but we do not check this here. '
        "We will now continue running with these files, as we can work with this here, "
        "but recommend to change the sample names.\n"
    )

del uniq_sample_names, problematic_filenames, relative_filenames
del numeric_sample_names
