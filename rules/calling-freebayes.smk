# =================================================================================================
#     Variant Calling
# =================================================================================================

def know_variants_extra():
    if config["data"]["reference"].get("known-variants"):
        return " --haplotype-basis-alleles " + config["data"]["reference"].get("known-variants")
    else:
        return ""

rule call_variants:
    input:
        # Our custom script needs the genome and its fai file.
        ref=config["data"]["reference"]["genome"],
        fai=get_fai(),

        # Get the bam and bai files for all files. If this is too slow, we need to split,
        # similar to what our GATK HaplotypeCaller rule does.
        # Without bai files, freebayes claims that it recomputes them, but actually crashes...
        samples=get_all_bams(),
        indices=get_all_bais(),

        # If known variants are set, use this as input to ensure the file exists.
        # We use the file via the know_variants_extra() function call below.
        known=config["data"]["reference"].get("known-variants"), # empty if key not present

        # If we restict the calling to some regions, use this file here.
        regions="called/{contig}.regions.bed" if config["settings"].get("restrict-regions") else []
    output:
        pipe("called/{contig}.vcf")
    log:
        "logs/freebayes/{contig}.log"
    benchmark:
        "benchmarks/freebayes/{contig}.bench.log"
    params:
        # Optional extra parameters.
        extra=config["params"]["freebayes"]["extra"] + know_variants_extra(),

        # Reference genome chunk size for parallelization (default: 100000)
        chunksize=config["params"]["freebayes"]["chunksize"]
    threads:
        # Need to exclude threads that we need for compression
        max(
            1,
            int(config["params"]["freebayes"]["threads"]) -
            int(config["params"]["freebayes"]["compress-threads"])
        )
    group:
        "call_variants"
    # wrapper:
    #     "0.55.1/bio/freebayes"
    conda:
        "../envs/freebayes.yaml"
    script:
        # We use our own version of the wrapper here, which improves cluster throughput.
        "../scripts/freebayes.py"

# Picard does not understand the bcf files that freebayes produces, so we have to take
# an unfortunate detour via vcf, and compress on-the-fly using a piped rule from above.
rule compress_vcf:
    input:
        "called/{contig}.vcf"
    output:
        protected("called/{contig}.vcf.gz")
    log:
        "logs/compress_vcf/{contig}.log"
    threads:
        config["params"]["freebayes"]["compress-threads"]
    group:
        "call_variants"
    conda:
        "../envs/tabix.yaml"
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
