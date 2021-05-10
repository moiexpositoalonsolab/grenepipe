# =================================================================================================
#     Trimming
# =================================================================================================

# The fastp wrapper is different from the other trimming wrappers in the the snakemake wrapper
# respository, because consistency is just not their strength... So, we have to provide an extra
# function to unpack the fastq file names into a list for us.
def unpack_fastp_files(wildcards):
    return list(get_fastq(wildcards).values())

rule trim_reads_se:
    input:
        sample=unpack_fastp_files
    output:
        trimmed="trimmed/{sample}-{unit}.fastq.gz",
        html="trimmed/{sample}-{unit}-se-fastp.html",
        json="trimmed/{sample}-{unit}-se-fastp.json"
    log:
        "logs/fastp/{sample}-{unit}.log"
    benchmark:
        "benchmarks/fastp/{sample}-{unit}.bench.log"
    params:
        extra=config["params"]["fastp"]["se"]
    threads:
        config["params"]["fastp"]["threads"]
    wrapper:
        "0.64.0/bio/fastp"

rule trim_reads_pe:
    input:
        sample=unpack_fastp_files
    output:
        trimmed=["trimmed/{sample}-{unit}.1.fastq.gz", "trimmed/{sample}-{unit}.2.fastq.gz"],
        html="trimmed/{sample}-{unit}-pe-fastp.html",
        json="trimmed/{sample}-{unit}-pe-fastp.json"
    log:
        "logs/fastp/{sample}-{unit}.log"
    benchmark:
        "benchmarks/fastp/{sample}-{unit}.bench.log"
    params:
        extra=config["params"]["fastp"]["pe"]
    threads:
        config["params"]["fastp"]["threads"]
    wrapper:
        "0.64.0/bio/fastp"

rule trim_reads_pe_merged:
    input:
        sample=unpack_fastp_files
    output:
        # Need to leave "trimmed" empty here, so that the wrapper works properly with merged,
        # so we use "merged" instead, and use it as an extra param.
        merged="trimmed/{sample}-{unit}-merged.fastq.gz",
        html="trimmed/{sample}-{unit}-pe-merged-fastp.html",
        json="trimmed/{sample}-{unit}-pe-merged-fastp.json"
    log:
        "logs/fastp/{sample}-{unit}.log"
    benchmark:
        "benchmarks/fastp/{sample}-{unit}.bench.log"
    params:
        extra=config["params"]["fastp"]["pe"] +
            " --merge --merged_out trimmed/{sample}-{unit}-merged.fastq.gz" +
            " --out1 trimmed/{sample}-{unit}-unmerged.pass-1.fastq.gz" +
            " --out2 trimmed/{sample}-{unit}-unmerged.pass-2.fastq.gz" +
            " --unpaired1 trimmed/{sample}-{unit}-unmerged.unpaired-1.fastq.gz" +
            " --unpaired2 trimmed/{sample}-{unit}-unmerged.unpaired-2.fastq.gz"
    threads:
        config["params"]["fastp"]["threads"]
    wrapper:
        "0.64.0/bio/fastp" # this runs fastp 0.20.0

# =================================================================================================
#     Trimming Results
# =================================================================================================

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(samples.loc[(sample, unit), "fq2"])

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if is_single_end(**wildcards):
        # single end sample
        return [ "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards) ]
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return [ "trimmed/{sample}-{unit}-merged.fastq.gz".format(**wildcards) ]
    else:
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{pair}.fastq.gz", pair=[1, 2], **wildcards)

def get_trimming_report(sample, unit):
    """Get the report needed for MultiQC."""
    if is_single_end(sample, unit):
        # single end sample
        return "trimmed/" + sample + "-" + unit + "-se-fastp.json"
    elif config["settings"]["merge-paired-end-reads"]:
        # merged paired-end samples
        return "trimmed/" + sample + "-" + unit + "-pe-merged-fastp.json"
    else:
        # paired-end sample
        return "trimmed/" + sample + "-" + unit + "-pe-fastp.json"
