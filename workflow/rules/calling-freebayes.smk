import platform

# =================================================================================================
#     Variant Calling
# =================================================================================================


def know_variants_extra():
    if config["data"]["known-variants"]:
        return " --haplotype-basis-alleles " + config["data"]["known-variants"]
    else:
        return ""


# Need to calculate the threads to be used for freebayes here already,
# see https://github.com/snakemake/snakefmt/issues/240
# Need to exclude threads that we need for compression:
freebayes_threads = max(
    1,
    int(config["params"]["freebayes"]["threads"])
    - int(config["params"]["freebayes"]["compress-threads"]),
)


rule call_variants:
    input:
        # Our custom script needs the genome, with index files, and its fai file.
        # The fai is integrated per checkpoint, which for this step here is not needed,
        # but is needed for the parallelization over configs, see the rule merge_variants.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        fai=get_fai,
        # Get the bam and bai files for all files. If this is too slow, we need to split,
        # similar to what our GATK HaplotypeCaller rule does.
        # Without bai files, freebayes claims that it recomputes them, but actually crashes...
        samples=get_all_bams(),
        indices=get_all_bais(),
        # If known variants are set, use this as input to ensure the file exists.
        # We use the file via the know_variants_extra() function call below,
        # but request it here as well to ensure that it and its index are present.
        known=config["data"]["known-variants"],
        knownidx=(
            config["data"]["known-variants"] + ".tbi" if config["data"]["known-variants"] else []
        ),
        # If we restict the calling to some regions, use this file here.
        # regions="calling/regions/{contig}.bed" if config["settings"].get("restrict-regions") else []
        regions=(
            "calling/regions/{contig}.bed"
            if (config["settings"].get("restrict-regions"))
            else (
                "calling/contig-groups/{contig}.bed"
                if (config["settings"].get("contig-group-size"))
                else []
            )
        ),
    output:
        pipe("calling/called/{contig}.vcf"),
        touch("calling/called/{contig}.vcf.done"),
    log:
        "logs/calling/freebayes/{contig}.log",
    benchmark:
        "benchmarks/calling/called/freebayes/{contig}.log"
    params:
        # Optional extra parameters.
        extra=config["params"]["freebayes"]["extra"] + know_variants_extra(),
        # Reference genome chunk size for parallelization (default: 100000)
        chunksize=config["params"]["freebayes"]["chunksize"],
    threads: freebayes_threads
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
        "calling/called/{contig}.vcf",
    output:
        (
            "calling/called/{contig}.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/called/{contig}.vcf.gz")
        ),
        # protected("calling/called/{contig}.vcf.gz")
        touch("calling/called/{contig}.vcf.gz.done"),
    log:
        "logs/calling/compress_vcf/{contig}.log",
    threads: config["params"]["freebayes"]["compress-threads"]
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


# Need an input function to work with the fai checkpoint
def merge_variants_vcfs_input(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("calling/called/{contig}.vcf.gz", contig=get_contigs(fai))


rule merge_variants:
    input:
        # fai is needed to calculate aggregation over contigs below.
        # This is the step where the genome is split into its contigs for parallel execution.
        # The get_fai() function uses a snakemake checkpoint to make sure that the fai is
        # produced before we use it here to get its content.
        ref=get_fai,
        contig_groups=contigs_groups_input,
        # The wrapper expects input to be called `vcfs`, but we can use `vcf.gz` as well.
        # vcfs=lambda w: expand("calling/called/{contig}.vcf.gz", contig=get_contigs())
        vcfs=merge_variants_vcfs_input,
    output:
        vcf="calling/genotyped-all.vcf.gz",
        done=touch("calling/genotyped-all.done"),
    params:
        # See duplicates-picard.smk for the reason whe need this on MacOS.
        extra=(
            " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true" if platform.system() == "Darwin" else ""
        ),
    log:
        "logs/calling/picard/merge-genotyped.log",
    benchmark:
        "benchmarks/calling/genotyped/picard/merge-genotyped.log"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "0.51.3/bio/picard/mergevcfs"
