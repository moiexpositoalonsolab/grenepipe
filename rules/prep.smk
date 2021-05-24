# =================================================================================================
#     Reference Genome
# =================================================================================================

# Helper function to get the name of the genome dictorary file as expected by GATK
def genome_dict():
    return os.path.splitext(config["data"]["reference"]["genome"])[0] + ".dict"

# We define some local variables for simplicity, and delete them later
# in order to not spam the global scope by accident.

# Get file names from config file. The reference genome file has already been stripped of the
# `.gz` extension if present in common.
genome=config["data"]["reference"]["genome"]

# We need to remove absolute paths here, otherwise the log files will contain broken paths.
genomename = os.path.basename(config["data"]["reference"]["genome"])

# In all rules below, we use hard coded file names (no wildcards), as snakemake cannot handle
# absolute file paths properly and gives us no reasonable way to use lambdas in the `log` part
# of the rule, which hence would lead to log file paths containing the absolute file path
# of our input genome. We do not want that - and as this whole prep script here only serves
# one purpose (prepare one genome for a given config file), we just hard code for simplicity.

# Write fai indices for the fasta reference genome file.
# This is a checkpoint, as downstream rules that parallelize over chromosomes/contigs of the
# reference genome will need to read this file in order to get the list of contigs.
# See get_fai() in calling.smk for the usage of this checkpoint.
checkpoint samtools_faidx:
    input:
        genome
    output:
        genome + ".fai"
    log:
        "logs/" + genomename + ".samtools_faidx.log"
    params:
        "" # optional params string
    wrapper:
        "0.51.3/bio/samtools/faidx"

# Uncompress the reference genome if it is gz, without deleting the original.
# Also, we provide our test data in gz-compressed form, in order to keep data in the git repo low.
# Hence, we have to decompress first.
rule decompress_genome:
    input:
        genome + ".gz"
    output:
        genome
    log:
        "logs/" + genomename + ".decompress.log"
    shell:
        "gunzip --keep {input}"

localrules:
    decompress_genome

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
        "logs/" + genomename + ".bwa_index.log"
    params:
        prefix=genome,
        algorithm="bwtsw"
    wrapper:
        "0.51.3/bio/bwa/index"

# Write a dictionary file for the genome.
# The input file extension is replaced by `dict`, instead of adding to it, so we have to trick
# around with the output file name here.
rule sequence_dictionary:
    input:
        genome
    output:
        genome_dict()
    # params:
    #     base= lambda wc: os.path.splitext(genome)[0],
    log:
        "logs/" + genomename + ".sequence_dictionary.log"
    conda:
        "../envs/prep.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output} > {log} 2>&1"

# Clean up the variables that we used above
del genome
del genomename

# =================================================================================================
#     Known Variants
# =================================================================================================

# We need the variants entry in the config to be either an empty list or a file path, which we
# already ensure in common.smk. This is becausae snakemake does not accept empty strings as input
# files - it has to either a file path or an empty list. So here we need to do a bit of trickery
# to allow for the case that no known variants file is given in the config:
# We define a local variable that is always a string. If it is empty, that does not seem to matter
# here, because those rules will never be invoked in that case, and so, snakemake does not seem to
# fail then. Still, we make it even more fail safe by setting it to a dummy string then that
# will just lead to rules that are never executed (in the case of no known variants file).
variants=config["data"]["reference"]["known-variants"]
if isinstance(variants, list) or not variants:
    if len(variants) > 0:
        raise Exception("Known variants has to be either a file path or an empty list." )
    variants="dummyfile"
else:
    # Somehow, some tool (was it GATK?) requires known variants to be in vcf.gz format,
    # so let's ensure this, and overwrite the config. We then also set our local variants to
    # the file _without_ the gz extension, so that we can set up the rules below correctly to
    # compress the file and build an index for it.
    if os.path.splitext(variants)[1] == ".vcf":
        # Set the config, so that the rules that actually use this file request the correct one.
        config["data"]["reference"]["known-variants"] += ".gz"
    elif variants.endswith(".vcf.gz"):
        # Set the local one to without the extension, to keep our rules below simple.
        variants = os.path.splitext(variants)[0]
    else:
        raise Exception(
            "Invalid known variants file type: '" + variants +
            "'. Needs to be either .vcf or .vcf.gz for some of the tools to work."
        )

# Compress the known variants vcf file using gzip, as this seems needed for GATK.
rule variants_vcf_compress:
    input:
        variants
    output:
        variants + ".gz"
    group:
        "known_variants"
    log:
        os.path.dirname(variants) + "/logs/" + os.path.basename(variants) + ".vcf_compress.log"
    wrapper:
        "0.27.1/bio/vcf/compress"

# Write an index file for the known variants file.
rule variants_vcf_index:
    input:
        variants + ".gz"
    output:
        variants + ".gz.tbi"
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf"
    group:
        "known_variants"
    log:
        os.path.dirname(variants) + "/logs/" + os.path.basename(variants) + ".vcf_index.log"
    wrapper:
        "0.55.1/bio/tabix"

# Clean up.
del variants
