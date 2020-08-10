# =================================================================================================
#     Variant Calling
# =================================================================================================

rule call_variants:
    input:
        ref=config["data"]["reference"]["genome"],

        # Get the bam and bai files for all files. If this is too slow, we need to split,
        # similar to what our GATK HaplotypeCaller rule does.
        # Without bai files, freebays claims that it recomputes them, but actually crashes...
        samples=get_all_bams(),
        indices=get_all_bais(),

        # If we use restricted regions, set them here
        regions="called/{contig}.regions.bed" if config["settings"].get("restrict-regions") else []
    output:
        pipe("called/{contig}.vcf")
    log:
        "logs/freebayes/{contig}.log"
    benchmark:
        "benchmarks/freebayes/{contig}.bench.log"
    params:
        extra=config["params"]["freebayes"]["extra"],         # optional parameters
        chunksize=config["params"]["freebayes"]["chunksize"]  # reference genome chunk size for parallelization (default: 100000)
    threads:
        config["params"]["freebayes"]["threads"]
    wrapper:
        "0.55.1/bio/freebayes"

# Picard does not understand the bcf files that freebayes produces, so we have to take
# an unfortunate detour via vcf, and compress on-the-fly using a piped rule from above.
rule compress_vcf:
    input:
        "called/{contig}.vcf"
    output:
        protected("called/{contig}.vcf.gz")
    log:
        "logs/compress_vcf/{contig}.log"
    conda:
        "../envs/tabix.yaml"
    threads:
        config["params"]["freebayes"]["compress-threads"]
    shell:
        # We need to "force" overwriting the file, because snakemake creates a named pipe,
        # which then already exists once bgzip gets active.
        "bgzip --force --threads {threads} {input[0]} > {output[0]} 2> {log}"

# =================================================================================================
#     Combining Calls
# =================================================================================================

rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below

        # The wrapper expects input to be called `vcfs`, but we can use `vcf.gz` as well.
        vcfs=lambda w: expand("called/{contig}.vcf.gz", contig=get_contigs())
    output:
        vcf="genotyped/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    benchmark:
        "benchmarks/picard/merge-genotyped.bench.log"
    wrapper:
        "0.51.3/bio/picard/mergevcfs"
