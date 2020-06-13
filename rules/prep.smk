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
    variants=config["data"]["reference"]["known-variants"]
    if os.path.splitext(variants)[1] == ".vcf":
        variants += ".gz"
    elif not variants.endswith(".vcf.gz"):
        raise Exception("Invalid known variants file type: " + variants )

genomedir=os.path.dirname(config["data"]["reference"]["genome"]) + "/"

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
        dict=os.path.splitext(genome)[0] + ".dict",
        vcf=variants
    group:
        "prep"

# All of the prep ruls are local. No need to submit 1min jobs to the cluster.
localrules: preparation, decompress_genome, bwa_index, samtools_faidx, sequence_dictionary, compress_vcf

# =================================================================================================
#     Prepare Genome
# =================================================================================================

# We provide our test data in gz-compressed form, in order to keep data in the git repo low.
# Hence, we have to decompress first.
rule decompress_genome:
    input:
        genome + ".gz"
    output:
        genome
    log:
        genomedir + "logs/" + genome + ".decompress.log"
    group:
        "prep"
    shell:
        "gunzip {input}"

# Write indices for a given fasta reference genome file
rule bwa_index:
    input:
        "{genome}"
    output:
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt",
        "{genome}.pac",
        "{genome}.sa"
    log:
        genomedir + "logs/{genome}.bwa_index.log"
    group:
        "prep"
    params:
        prefix="{genome}",
        algorithm="bwtsw"
    wrapper:
        "0.51.3/bio/bwa/index"

# Write more indices for the fasta reference genome file
rule samtools_faidx:
    input:
        "{genome}"
    output:
        "{genome}.fai"
    log:
        genomedir + "logs/{genome}.samtools_faidx.log"
    group:
        "prep"
    params:
        "" # optional params string
    wrapper:
        "0.51.3/bio/samtools/faidx"

# Write a dictionary file for the genome.
# This file is expected to replace the file extension, instead of adding to it.
# We also add the genome file itself as an input, so that it is decompressed.
rule sequence_dictionary:
    input:
        base="{genome}" + os.path.splitext(genome)[1],
        file=genome
    output:
        "{genome}.dict"
    log:
        genomedir + "logs/{genome}.sequence_dictionary.log"
    group:
        "prep"
    conda:
        "../envs/prep.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input.base} -O {output} > {log} 2>&1"

# Compress a vcf file using gzip
rule compress_vcf:
    input:
        "{prefix}.vcf"
    output:
        "{prefix}.vcf.gz"
    log:
        genomedir + "logs/{prefix}.compress_vcf.log"
    group:
        "prep"
    wrapper:
        "0.27.1/bio/vcf/compress"
