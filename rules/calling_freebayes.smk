# =================================================================================================
#     Helper Functions
# =================================================================================================

from itertools import chain

# Return the bam file(s) for a given sample
def get_sample_bams(sample):
    return expand(get_mapping_result(), sample=sample, unit=samples.loc[sample].unit)

# Return the bai file(s) for a given sample
def get_sample_bais(sample):
    return expand(get_mapping_result(True), sample=sample, unit=samples.loc[sample].unit)

# Return the bam file(s) for all samples
def get_all_bams():
    return list(chain.from_iterable( [ get_sample_bams(sample) for sample in sample_names ] ))

# Return the bai file(s) for all samples
def get_all_bais():
    return list(chain.from_iterable( [ get_sample_bais(sample) for sample in sample_names ] ))

# =================================================================================================
#     Variant Calling
# =================================================================================================

rule call_variants:
    input:
        ref=config["data"]["reference"]["genome"],

        # Get the bam and bai files for the given sample only.
        # Without bai files, freebays claims that it recomputes them, but actually crashes...
        samples=get_all_bams(),
        indices=get_all_bais(),

        # If we use restricted regions, set them here
        regions=config["rundir"] + "called/{contig}.regions.bed" if config["settings"].get("restrict-regions") else []
    output:
        pipe(config["rundir"] + "called/{contig}.vcf")
    log:
        config["rundir"] + "logs/freebayes/{contig}.log"
    params:
        extra=config["params"]["freebayes"]["extra"],         # optional parameters
        chunksize=config["params"]["freebayes"]["chunksize"]  # reference genome chunk size for parallelization (default: 100000)
    threads:
        4 # TODO
    wrapper:
        "0.55.1/bio/freebayes"

# Picard does not understand the bcf files that freebayes produces, so we have to take
# an unfortunate detour via vcf, and compress on-the-fly using a piped rule.
rule compress_vcf:
    input:
        config["rundir"] + "called/{contig}.vcf"
    output:
        protected(config["rundir"] + "called/{contig}.vcf.gz")
    log:
        config["rundir"] + "logs/compress_vcf/{contig}.log"
    conda:
        "../envs/tabix.yaml"
    threads:
        2 # TODO
    shell:
        # We need to "force" overwriting the file, because snakemake creates a named pipe,
        # which then already exists once bgzip gets active.
        "bgzip --force --threads {threads} {input[0]} > {output[0]} 2> {log} "

# =================================================================================================
#     Combining Calls
# =================================================================================================

rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below

        # The wrapper expects input to be called `vcfs`, but we can use `vcf.gz` as well.
        vcfs=lambda w: expand(config["rundir"] + "called/{contig}.vcf.gz", contig=get_contigs())
    output:
        vcf=config["rundir"] + "genotyped/all.vcf.gz"
    log:
        config["rundir"] + "logs/picard/merge-genotyped.log"
    wrapper:
        "0.51.3/bio/picard/mergevcfs"
