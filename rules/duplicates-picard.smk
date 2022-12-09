import platform

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
        metrics="qc/dedup/{sample}-{unit}.metrics.txt",
        done=touch("dedup/{sample}-{unit}.done")
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
    benchmark:
        "benchmarks/picard/dedup/{sample}-{unit}.bench.log"
    params:
        # Take the params from the config.
        # On MacOS (we experienced it with 10.16, 11, and 12 so far), there is an issue between Java
        # and some system libraries, leading the JRE to exit with fatal error SIGSEGV caused by
        # libgkl_compression, see https://github.com/broadinstitute/picard/issues/1329.
        # Hence, on MacOS, we add the two settings recommended by the github issue.
        config["params"]["picard"]["MarkDuplicates"] + (
            " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"
            if platform.system() == "Darwin"
            else ""
        )
    group:
        "mapping_extra"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "0.51.3/bio/picard/markduplicates"
