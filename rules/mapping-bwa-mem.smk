# =================================================================================================
#     Read Mapping
# =================================================================================================

rule map_reads:
    input:
        reads=get_trimmed_reads
    output:
        config["rundir"] + "mapped/{sample}-{unit}.sorted.bam"
    log:
        config["rundir"] + "logs/bwa-mem/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/bwa-mem/{sample}-{unit}.bench.log"
    params:
        index=config["data"]["reference"]["genome"],

        # We need the read group tags, including `ID` and `SM`, as downstream tools use these.
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        # TODO add LD field as well for the unit?! http://www.htslib.org/workflow/

        sort="samtools",         # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate", # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads:
        config["params"]["bwamem"]["threads"]
    resources:
        # Increase time limit in factors of 2h, if the job fails due to time limit.
        time = lambda wildcards, attempt: int(120 * int(attempt))
    wrapper:
        "0.51.3/bio/bwa/mem"
