# =================================================================================================
#     Dummy Trimming
# =================================================================================================


def get_trimmed_reads(wildcards):
    # Simply forward the files to their original fastq files.
    return list(get_fastq(wildcards).values())


# In our current setup, we require "done" files to trick snakemake into working properly.
# These might not exist when running with raw fastq files, so we touch them here.
def get_trimmed_reads_done(wildcards):
    files = get_trimmed_reads(wildcards)

    # Touch all non-existing files. If they already exist,
    # we do nothing, to not mess with their time stamps.
    for f in files:
        if not os.path.isfile(f):
            Path(f + ".done").touch()

    # Now we can return the fastq done file list to the caller.
    return [f + ".done" for f in files]


def get_trimming_report(sample, unit):
    # No trimming report here.
    return []
