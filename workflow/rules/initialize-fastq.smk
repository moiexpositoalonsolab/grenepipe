# =================================================================================================
#     Read Samples Table
# =================================================================================================

# We add the samples information to the config, in order to not spam our global scope
# (unfortunately, while python is super good with namespaces, it is super bad with scopes,
# and in particular in snakemake, every file is included so that all variables defined in a file
# are global and accessible in all subsequent files as well...).
# We use a new top level key that is not used in the config file for this. Assert this.
if "global" in config:
    raise Exception("Config key 'global' already defined. Someone messed with our setup.")
else:
    config["global"] = {}

# Read samples and units table, and enforce to use strings in the index
config["global"]["samples"] = pd.read_csv(
    config["data"]["samples-table"], sep="\t", dtype=str
).set_index(["sample", "unit"], drop=False)
config["global"]["samples"].index = config["global"]["samples"].index.set_levels(
    [i.astype(str) for i in config["global"]["samples"].index.levels]
)
snakemake.utils.validate(config["global"]["samples"], schema="../schemas/samples.schema.yaml")


# Get a list of all samples names, in the same order as the input sample table.
# Samples with multiple units appear only once, at the first position in the table.
# We cannot use a simple approach here, as this messes up the sample
# order, which we do not want... (good that we noticed that bug though!)
# So instead, we iterate, and add sample names incrementally.
config["global"]["sample-names"] = list()
for index, row in config["global"]["samples"].iterrows():
    s = row["sample"]
    if s not in config["global"]["sample-names"]:
        config["global"]["sample-names"].append(s)

# Unordered list of all unit names that appear in all the samples.
config["global"]["unit-names"] = list(
    set(config["global"]["samples"].index.get_level_values("unit"))
)


# Helper function to get a list of all units of a given sample name.
def get_sample_units(sample):
    res = list()
    for unit in config["global"]["samples"].loc[sample].unit:
        if unit not in res:
            res.append(unit)
    return res


# Helper function for user output to summarize the samples found in the table.
def get_sample_units_print():
    unitcnt = len(config["global"]["samples"].index.get_level_values("unit"))
    if unitcnt == len(config["global"]["sample-names"]):
        return str(len(config["global"]["sample-names"]))
    else:
        return (
            str(len(config["global"]["sample-names"])) + ", with " + str(unitcnt) + " total units"
        )


# Wildcard constraints: only allow sample names from the spreadsheet to be used
wildcard_constraints:
    sample="|".join(config["global"]["sample-names"]),
    unit="|".join(config["global"]["unit-names"]),


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
if not os.path.isabs(config["data"]["samples-table"]):
    logger.warning(
        "Path to the samples table as provided in the config file is not an absolute path. "
        "We recommend using absolute paths for all files.\n"
    )

# List that contains tuples for all samples with their units.
# In other words, a list of tuples of the sample and unit column of the sample table,
# in the same order.
config["global"]["sample-units"] = list()
problematic_filenames = 0
relative_filenames = 0
for index, row in config["global"]["samples"].iterrows():
    if (row["sample"], row["unit"]) in config["global"]["sample-units"]:
        raise Exception(
            "Multiple rows with identical sample name and unit found in samples table: "
            + str(row["sample"])
            + " "
            + str(row["unit"])
        )
    config["global"]["sample-units"].append((row["sample"], row["unit"]))

    # Do a check that the sample and unit names are valid file names.
    # They are used for file names, and would cause weird errors if they contain bad chars.
    if not valid_filename(row["sample"]) or not valid_filename(row["unit"]):
        raise Exception(
            "Invalid sample name or unit name found in samples table that contains characters "
            + "which cannot be used as sample/unit names for naming output files: "
            + str(row["sample"])
            + " "
            + str(row["unit"])
            + "; for maximum robustness, we only allow alpha-numerical, dots, dashes, "
            "and underscores. "
            + "Use for example the script `tools/copy-samples.py --mode link [...] --clean` "
            + "to create a new samples table and symlinks to the existing fastq files to solve this."
        )

    # Do a check of the fastq file names.
    if not os.path.isfile(row["fq1"]) or (not pd.isnull(row["fq2"]) and not os.path.isfile(row["fq2"])):
        raise Exception(
            "Input fastq files listed in the input files samples table "
            + config["data"]["samples-table"]
            + " not found: "
            + str(row["fq1"])
            + "; "
            + str(row["fq2"])
        )
    if not valid_filepath(row["fq1"]) or (not pd.isnull(row["fq2"]) and not valid_filepath(row["fq2"])):
        problematic_filenames += 1
    if not os.path.isabs(row["fq1"]) or (not pd.isnull(row["fq2"]) and not os.path.isabs(row["fq2"])):
        relative_filenames += 1

# Warning about input names and files.
if problematic_filenames > 0:
    logger.warning(
        str(problematic_filenames)
        + " of the "
        + str(len(config["global"]["sample-names"]))
        + " input samples listed in the input files samples table "
        + config["data"]["samples-table"]
        + " contain problematic characters. We generally advise to only use alpha-numeric "
        "characters, dots, dashes, and underscores. "
        "Use for example the script `tools/copy-samples.py --mode link [...] --clean` "
        "to create a new samples table and symlinks to the existing fastq files to solve this. "
        "We will try to continue running with these files, but it might lead to errors.\n"
    )
if relative_filenames > 0:
    logger.warning(
        str(relative_filenames)
        + " of the "
        + str(len(config["global"]["sample-names"]))
        + " input samples listed in the input files samples table "
        + config["data"]["samples-table"]
        + " use relative file paths. We generally advise to only use absolute paths. "
        "Use for example the script `tools/generate-table.py` "
        "to create a samples table with absolute paths. "
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
        + " input sample names listed in the input files samples table "
        + config["data"]["samples-table"]
        + " are just numbers, or can be converted to numbers. We generally advise to avoid this, "
        "as it might confuse downstream processing such as Excel and R. The same applies for names "
        'that can be converted to dates ("Dec21"), but we do not check this here. '
        "We will now continue running with these files, as we can work with this here, "
        "but recommend to change the sample names.\n"
    )

del problematic_filenames, relative_filenames
del numeric_sample_names
