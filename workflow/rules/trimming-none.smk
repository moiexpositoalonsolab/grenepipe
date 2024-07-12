# =================================================================================================
#     Dummy Trimming
# =================================================================================================


def get_trimmed_reads(wildcards):
    # Simply forward the files to their original fastq files.
    return get_fastq(wildcards).values()


def get_trimming_report(sample, unit):
    # No trimming report here.
    return []
