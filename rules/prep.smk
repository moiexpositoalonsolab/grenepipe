# =================================================================================================
#     Reference Genome
# =================================================================================================

# Helper function to get the name of the genome dictorary file as expected by GATK
def genome_dict():
    return os.path.splitext(config["data"]["reference-genome"])[0] + ".dict"

# We define some local variables for simplicity, and delete them later
# in order to not spam the global scope by accident.

# Get file names from config file. The reference genome file has already been stripped of the
# `.gz` extension if present in common.
genome=config["data"]["reference-genome"]

# We need to remove absolute paths here, otherwise the log files will contain broken paths.
genomename = os.path.basename( config["data"]["reference-genome"] )
genomedir  = os.path.dirname(  config["data"]["reference-genome"] )

# We write the log file to where the known variants file is, so that this is independent
# of the particular run, making the log files easier to find for users.
genome_logdir = os.path.join( genomedir, "logs" )

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
        os.path.join( genome_logdir, genomename + ".samtools_faidx.log" )
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
        os.path.join( genome_logdir, genomename + ".decompress.log" )
    shell:
        # Cannot use gunzip here, as CentOS does not support the --keep option...
        # "gunzip --keep {input}"
        # Well, zcat doesn't work on MacOS, see https://serverfault.com/a/704521
        # "zcat {input} > {output}"
        # So back to gunzip, but with a different arg to keep the original file...
        "gunzip -c {input} > {output}"

localrules:
    decompress_genome

# Write indices for a given fasta reference genome file.
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
        os.path.join( genome_logdir, genomename + ".bwa_index.log" )
    params:
        prefix=genome,
        algorithm="bwtsw"
    wrapper:
        "0.51.3/bio/bwa/index"

# By default, we use the above bwa index mapping, but here we also have a rule for bwa mem2 indexing,
# as this needs some special indices. However, it unfortunately does not produce all indices
# of the original bwa index (because... bioinformatics...), so we have to work around this,
# without creating naming conflicts. Hence, when using bwa mem2, we create its extra indices
# in a temp dir, and then only move the ones that are not also produced by bwa index... so hacky.
if config["settings"]["mapping-tool"] == "bwamem2":

    # We do not use the wrapper here, as it would otherwise overwrite the indices created by
    # bwa index above... and working around that while using the wrapper is super hacky,
    # so instead we copy over the wrapper code and adapt it to our needs. The most important
    # adaptation is to sym-link to the ref genome in the temp dir, so that the extra indices
    # are created there... wow, hacky hacky.
    rule bwa_mem2_index:
        input:
            genome
        output:
            genome + ".0123",
            genome + ".bwt.2bit.64"

            # Files that are also created by bwa mem2 index: amb, ann, pac.
            # We checked with our test data, and they are identical to the ones produced by
            # bwa index above, so we can just delete them.
        params:
            tempdir  = genomedir + "/bwa-mem2-index-temp/",
            basename = genomedir + "/bwa-mem2-index-temp/" + genomename
        log:
            os.path.join( genome_logdir, genomename + ".bwa-mem2_index.log" )
        conda:
            "../envs/bwa-mem2.yaml"
        shell:
            # Somehow, we need all variables in params here,
            # as otherwise the commands below fail for some weird reason...
            "mkdir -p {params.tempdir} ; "
            "ln -sf {input} {params.basename} ; "
            "bwa-mem2 index {params.basename} > {log} 2>&1 ; "
            "mv {params.basename}.0123        {output[0]} ; "
            "mv {params.basename}.bwt.2bit.64 {output[1]} ; "
            "rm -r {params.tempdir}"

# Write a dictionary file for the genome.
# The input file extension is replaced by `dict`, instead of adding to it, so we have to trick
# around with the output file name here.
rule sequence_dictionary:
    input:
        genome
    output:
        genome_dict()
    params:
        # See duplicates-picard.smk for the reason whe need this on MacOS.
        extra = (
            " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
            if platform.system() == "Darwin"
            else ""
        )
    #     base= lambda wc: os.path.splitext(genome)[0],
    log:
        os.path.join( genome_logdir, genomename + ".sequence_dictionary.log" )
    conda:
        "../envs/picard.yaml"
    # We used GATK here for a while, but for whatever reason, in our GitHub Actions tests,
    # this just randomly failed to produce a valid dict file, with no indication as to why,
    # so we had to re-start CI jobs until it randomly worked... Trying Picard now, hoping to fix it.
    shell:
        "picard CreateSequenceDictionary REFERENCE={input} OUTPUT={output} {params.extra} > {log} 2>&1"
    # shell:
    #     "gatk CreateSequenceDictionary -R {input} -O {output} --VERBOSITY DEBUG > {log} 2>&1"

# Get some statistics about the reference genome
rule reference_seqkit:
    input:
        genome
    output:
        genome + ".seqkit"
    params:
        extra = config["params"]["seqkit"]["extra"]
    log:
        os.path.join( genome_logdir, genomename + ".seqkit.log" )
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit stats {input} {params.extra} > {output} 2> {log}"

# =================================================================================================
#     Known Variants
# =================================================================================================

# We need the variants entry in the config to be either an empty list or a file path, which we
# already ensure in common.smk. This is because snakemake does not accept empty strings as input
# files - it has to be either a file path or an empty list. So here we need to do a bit of trickery
# to allow for the case that no known variants file is given in the config:
# We define a local variable that is always a string. If it is empty, that does not seem to matter
# here, because those rules will never be invoked in that case, and so, snakemake does not seem to
# fail then. Still, we make it even more fail safe by setting it to a dummy string then that
# will just lead to rules that are never executed (in the case of no known variants file).
variants=config["data"]["known-variants"]
has_known_variants = True
if isinstance(variants, list) or not variants:
    if len(variants) > 0:
        raise Exception("Known variants has to be either a file path or an empty list." )
    variants="dummyfile"
    has_known_variants = False
else:
    # Somehow, some tool (was it GATK?) requires known variants to be in vcf.gz format,
    # so let's ensure this, and overwrite the config. We then also set our local variants to
    # the file _without_ the gz extension, so that we can set up the rules below correctly to
    # compress the file and build an index for it.
    if os.path.splitext(variants)[1] == ".vcf":
        # Set the config, so that the rules that actually use this file request the correct one.
        config["data"]["known-variants"] += ".gz"
    elif variants.endswith(".vcf.gz"):
        # Set the local one to without the extension, to keep our rules below simple.
        variants = os.path.splitext(variants)[0]
    else:
        raise Exception(
            "Invalid known variants file type: '" + variants +
            "'. Needs to be either .vcf or .vcf.gz for some of the tools to work."
        )

# We write the log file to where the known variants file is, so that this is independent
# of the particular run, making the log files easier to find for users.
variant_logdir = os.path.join( os.path.dirname(variants), "logs" )

# Compress the known variants vcf file using gzip, as this seems needed for GATK.
rule variants_vcf_compress:
    input:
        variants
    output:
        variants + ".gz"
    group:
        "known_variants"
    log:
        os.path.join( variant_logdir, os.path.basename(variants) + ".vcf_compress.log" )
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
        os.path.join( variant_logdir, os.path.basename(variants) + ".vcf_index.log" )
    conda:
        "../envs/tabix.yaml"
    wrapper:
        "0.55.1/bio/tabix"

# =================================================================================================
#     All prep rule
# =================================================================================================

# This alternative target rule executes all prep rules, so that we can get all indices of the
# reference etc. This is useful in settings where we have an unprepared reference, but want
# to run the pipeline multiple times already. If we did not prep the genome first, each pipeline
# would try to do that on its own, which might lead to conflicts and broken files.
rule all_prep:
    input:
        ref=genome,
        ref_idcs=expand(
            genome + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),
        ref_idcs2=expand(
            genome + ".{ext}",
            ext=[ "0123", "bwt.2bit.64" ]
        ) if config["settings"]["mapping-tool"] == "bwamem2" else [],
        ref_dict=genome_dict(),
        ref_stat=genome + ".seqkit",
        known_vars=variants + ".gz.tbi" if has_known_variants else []

localrules: all_prep

# Clean up the variables that we used above
del genome
del genomename
del genomedir
del variants
