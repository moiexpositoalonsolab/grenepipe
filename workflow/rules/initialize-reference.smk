# =================================================================================================
#     Reference Validity Checks
# =================================================================================================


# Helper to check that a string contains no invalid chars for file names.
# This is just the file, not its path! Slashes are considered invalid by this function.
def valid_filename(fn):
    # Only accept alnum, underscore, dash, and dot.
    return fn.replace("_", "").replace("-", "").replace(".", "").isalnum() and fn.isascii()
    # Bit more loose: allow all except Windows forbidden chars.
    # return not( True in [c in fn for c in "<>:\"/\\|?*"])


# Helper to check if a file path contains weird characters.
# We just want to warn about this for input fastq files, but still try to continue.
def valid_filepath(fn):
    # Only accept alnum, underscore, dash, dot, and slashes.
    clean = fn.replace("_", "").replace("-", "").replace(".", "").replace("/", "").replace("\\", "")
    return clean.isalnum() and clean.isascii()


# We need to clean up the file name for the reference genome. We normalize the file name to be
# without the gz ending. The prep rule decompress_genome provides the unzipped genome as needed.
if config["data"]["reference-genome"].endswith(".gz"):
    config["data"]["reference-genome"] = os.path.splitext(config["data"]["reference-genome"])[0]

# We recommend to use absolute paths. Check that for the reference genome.
if not os.path.isabs(config["data"]["reference-genome"]):
    logger.warning(
        "Path to the reference genome as provided in the config file is not an absolute path. "
        "We recommend using absolute paths for all files.\n"
    )

# GATK only accepts reference genomes if their file ending is `fa` or `fasta`, at least in the
# version that we are using. This has been updated later to also include `fas`, see here:
# https://github.com/samtools/htsjdk/commit/657b0a6076d84582a19b741fc28cd3c9a12384bf#diff-77d63fdf4a920459b3a44ead94ad979b93115eba0749baa10a694d061b9d6c1f
# It would of course be better to check the file contents instead of the extension, but okay...
# So here, we check this, in order to provide some better error messages for users,
# instead of having them run into cryptic log messages "File is not a supported reference file type".
# Add `".fas", ".fas.gz"` later if we upgrade GATK.
fastaexts = (".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fas", ".fas.gz")
if not config["data"]["reference-genome"].endswith(fastaexts):
    raise Exception(
        "Reference genome file path does not end in "
        + str(fastaexts)
        + ", which unfortunately is needed by GATK to be able to find the file. "
        + "Please rename the file and change the path in the config.yaml"
    )

# Some settings in the config file need to be converted from empty string to empty list, it seems,
# so that rules that use the files specified in these settings are working properly.
# Maybe there is some better way, but for now, this is working.
# We need to do this after the scheme validation, as we are introducing changes here
# that are not according to the scheme.
if "known-variants" not in config["data"] or not config["data"]["known-variants"]:
    config["data"]["known-variants"] = []
if "restrict-regions" not in config["settings"] or not config["settings"]["restrict-regions"]:
    config["settings"]["restrict-regions"] = []

# We recommend to use absolute paths. Check that for the known variants.
if config["data"]["known-variants"] and not os.path.isabs(config["data"]["known-variants"]):
    logger.warning(
        "Path to samples table as provided in the config file is not an absolute path. "
        "We recommend using absolute paths for all files.\n"
    )
