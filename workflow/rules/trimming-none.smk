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

    # Check if we have touched the fastq done files already
    if not hasattr(get_trimmed_reads_done, "done"):
        get_trimmed_reads_done.done = False

    # If not, touch all files, then set the internal flag
    # so that we do not do this every time this function is called.
    if not get_trimmed_reads_done.done:
        for f in files:
            Path(f + ".done").touch()
    get_trimmed_reads_done.done = True

    # Now we can return the fastq done file list to the caller.
    return [f + ".done" for f in files]


def get_trimming_report(sample, unit):
    # No trimming report here.
    return []
