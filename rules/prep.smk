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

# Get file names from config file. The reference genome file has already been stripped of the
# `.gz` extension if present in common.
genome=config["data"]["reference"]["genome"]
if not config["data"]["reference"]["known-variants"]:
    variants=""
    variants_index=""
else:
    variants=config["data"]["reference"]["known-variants"]
    if os.path.splitext(variants)[1] == ".vcf":
        variants += ".gz"
    elif not variants.endswith(".vcf.gz"):
        raise Exception("Invalid known variants file type: " + variants )
    variants_index = variants + ".tbi"

# We need to remove absolute paths here, otherwise the log files will contain broken paths.
genomedir  = os.path.dirname(config["data"]["reference"]["genome"]) + "/"
genomename = os.path.basename(config["data"]["reference"]["genome"])

# The snake strikes again. For the "all" preparation rule below, we need `variants` to be either
# an actual file path, or an empty list, as Snakemake does not recognize empty strings properly...
# However, for the rules that work on the variants, we cannot use a list, so we need to provide
# (empty) strings in these cases. So ugly.
variantsdir  = os.path.dirname(variants) + "/"
variantsname = os.path.basename(variants)
if variants:
    variants_target=variants
    variants_index_target=variants_index
else:
    variants_target=[]
    variants_index_target=[]

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
        vcf=variants_target,
        vcfi=variants_index_target
    group:
        "prep"

# All of the prep ruls are local. No need to submit 1min jobs to the cluster.
localrules: preparation, decompress_genome, bwa_index, samtools_faidx, sequence_dictionary, vcf_compress, vcf_index

# =================================================================================================
#     Prepare Genome
# =================================================================================================

# In all rules below, we use hard coded file names (no wildcards), as stupid snakemake cannot
# handle absolute file paths properly and gives us no reasonable way to use lambdas in the `log`
# part of the rule, which hence would lead to log file paths containing the absolute file path
# of our input genome. We do not want that - and as this whole prep script here only serves
# one purpose (prepare one genome for a given config file), we just hard code for simplicity.

# We provide our test data in gz-compressed form, in order to keep data in the git repo low.
# Hence, we have to decompress first.
rule decompress_genome:
    input:
        genome + ".gz"
    output:
        genome
    log:
        genomedir + "logs/" + genomename + ".decompress.log"
    group:
        "prep"
    shell:
        "gunzip {input}"

# Write indices for a given fasta reference genome file
rule bwa_index:
    input:
        genome
    output:
        genome + ".amb",
        genome + ".ann",
        genome + ".bwt",
        genome + ".pac",
        genome + ".sa"
    log:
        genomedir + "logs/" + genomename + ".bwa_index.log"
    group:
        "prep"
    params:
        prefix=genome,
        algorithm="bwtsw"
    wrapper:
        "0.51.3/bio/bwa/index"

# Write more indices for the fasta reference genome file
rule samtools_faidx:
    input:
        genome
    output:
        genome + ".fai"
    log:
        genomedir + "logs/" + genomename + ".samtools_faidx.log"
    group:
        "prep"
    params:
        "" # optional params string
    wrapper:
        "0.51.3/bio/samtools/faidx"

# Write a dictionary file for the genome.
# The input file extension is replaced by `dict`, instead of adding to it, so we have to trick
# around with the output file name here.
rule sequence_dictionary:
    input:
        genome
    output:
        os.path.splitext(genome)[0] + ".dict"
    # params:
    #     base= lambda wc: os.path.splitext(genome)[0],
    log:
        genomedir + "logs/" + genomename + ".sequence_dictionary.log"
    group:
        "prep"
    conda:
        "../envs/prep.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output} > {log} 2>&1"

# Compress a vcf file using gzip
rule vcf_compress:
    input:
        variants
    output:
        variants + ".gz"
    log:
        genomedir + "logs/" + variantsname + ".vcf_compress.log"
    group:
        "prep"
    wrapper:
        "0.27.1/bio/vcf/compress"

rule vcf_index:
    input:
        variants
    output:
        variants + ".tbi"
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf"
    log:
        genomedir + "logs/" + variantsname + ".vcf_index.log"
    wrapper:
        "0.55.1/bio/tabix"
