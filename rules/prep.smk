import os

# Contains steps that are not part of the original Snakemake demo pipeline, but we incldue them here
# for full comfort. For the original steps, see
# https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/blob/master/.github/workflows/main.yml

# Unfortunately, at the moment, this Snakefile is not part of the main pipeline,
# as the DAG of the main pipeline depends on the output of this script.
# This could be solved with snakemake checkpoints, but for now, it is easier to run it separately:
# From the main directory: `snakemake --snakefile rules/prep.smk`

# =================================================================================================
#     Setup
# =================================================================================================

# We need to load the config file
include: "common.smk"

# Get file names from config file
genome=config["data"]["reference"]["genome"]
if not config["data"]["reference"]["known-variants"]:
    variants=[]
else:
    variants=config["data"]["reference"]["known-variants"] + ".gz"

# =================================================================================================
#     Main Rule
# =================================================================================================

rule preparation:
    input:
        amb=genome + ".amb",
        ann=genome + ".ann",
        bwt=genome + ".bwt",
        pac=genome + ".pac",
        sa=genome + ".sa",
        fai=genome + ".fai",
        # fadict=genome + ".dict",
        dict=os.path.splitext(genome)[0] + ".dict",
        vcf=variants

# =================================================================================================
#     Rules for preparing input files
# =================================================================================================

rule bwa_index:
    input:
        "{genome}"
        # genome=config["data"]["reference"]["genome"]
    output:
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt",
        "{genome}.pac",
        "{genome}.sa"
    log:
        "logs/bwa_index/{genome}.log"
    params:
        prefix="{genome}",
        algorithm="bwtsw"
    wrapper:
        "0.51.3/bio/bwa/index"

rule samtools_faidx:
    input:
        "{genome}"
    output:
        "{genome}.fai"
    params:
        "" # optional params string
    wrapper:
        "0.51.3/bio/samtools/faidx"

rule sequence_dictionary:
    input:
        "{genome}" + os.path.splitext(genome)[1]
    output:
        "{genome}.dict"
    conda:
        "../envs/prep.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"

rule compress_vcf:
    input:
        "{prefix}.vcf"
    output:
        "{prefix}.vcf.gz"
    wrapper:
        "0.27.1/bio/vcf/compress"
