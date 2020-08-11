# =================================================================================================
#     Mark Duplicates
# =================================================================================================

# Used in `mapping.smk`

rule mark_duplicates:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam="dedup/{sample}-{unit}.bam",
        metrics="qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
    benchmark:
        "benchmarks/picard/dedup/{sample}-{unit}.bench.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    group:
        "mapping-extra"
    wrapper:
        "0.51.3/bio/picard/markduplicates"
