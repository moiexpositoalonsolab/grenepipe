# =================================================================================================
#     Mark Duplicates
# =================================================================================================

# Used in `mapping.smk`

rule mark_duplicates:
    input:
        # Get either the normal mapping output, or, if additional filtering via `samtools view`
        # is set in the config settings: filter-mapped-reads, use the filtered output instead.
        get_mapped_reads
        # "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=(
            "dedup/{sample}-{unit}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("dedup/{sample}-{unit}.bam")
        ),
        metrics="qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
    benchmark:
        "benchmarks/picard/dedup/{sample}-{unit}.bench.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    group:
        "mapping_extra"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "0.51.3/bio/picard/markduplicates"
