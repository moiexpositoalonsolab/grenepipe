import os

# =================================================================================================
#     HAFpipe Task 1:  Make SNP Table
# =================================================================================================

# We get the list of chromosomes from the config that the user wants HAFpipe to run for,
# and cross check them with the actual ref genome fai file, to avoid irritation.
# We use wildcard {chrom} here on purpose, instead of {contig} that we use for the rest of grenepipe,
# in order to (a) make it clear that these need to be actual sequences from the reference genome,
# and not our contig groups for example, and to (b) avoid accidents when matching the wildcards.
def get_hafpipe_chromosomes( fai ):
    ref_chrs = get_chromosomes( fai )
    haf_chrs = list( config["params"]["hafpipe"]["chromosomes"] )
    haf_chrs = [ str(v) for v in haf_chrs ]
    for chr in haf_chrs:
        if not chr in ref_chrs:
            raise Exception(
                "Chromosome '" + chr + "' specified via the config `params: hafpipe: chromosomes` " +
                "list for running HAFpipe is not part of the reference genome."
            )
    return haf_chrs

# We allow users to specify a directory for the snp table, to avoid recomputation of Task 1.
def get_hafpipe_snp_table_dir():
    cfg_dir = config["params"]["hafpipe"]["snp-table-dir"]
    if cfg_dir:
        return cfg_dir.rstrip('/')
    return "hafpipe/snp-tables"

rule hafpipe_snp_table:
    input:
        vcf=config["params"]["hafpipe"]["founder-vcf"]
    output:
        snptable=get_hafpipe_snp_table_dir() + "/{chrom}.csv"
    params:
        tasks="1",
        chrom="{chrom}",
        extra=config["params"]["hafpipe"]["snp-table-extra"]
    log:
        "logs/hafpipe/snp-table/{chrom}.log"
    conda:
        "../envs/hafpipe.yaml"
    script:
        "../scripts/hafpipe.py"

# Get the list of snp table files.
def get_all_hafpipe_raw_snp_tables(wildcards):
    # We use the fai file of the ref genome to cross-check the list of chromosomes in the config.
    # We use a checkpoint to create the fai file from our ref genome, which gives us the chrom names.
    # Snakemake then needs an input function to work with the fai checkpoint here.
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand(
        get_hafpipe_snp_table_dir() + "/{chrom}.csv",
        chrom=get_hafpipe_chromosomes( fai )
    )

# Rule that requests all HAFpipe SNP table files, so that users can impute them themselves.
rule all_hafpipe_snp_tables:
    input:
        get_all_hafpipe_raw_snp_tables

localrules: all_hafpipe_snp_tables

# =================================================================================================
#     HAFpipe Task 2:  Impute SNP Table
# =================================================================================================

# Edge case: HAFpipe allows "none" instead of empty for the impmethod.
# We get rid of this here for simplicity.
if config["params"]["hafpipe"]["impmethod"] == "none":
    config["params"]["hafpipe"]["impmethod"] = ""

# Shorthand
impmethod = config["params"]["hafpipe"]["impmethod"]

# We want to distinguish between the two impute methods that HAFpipe offers,
# and the case that the user provided their own script via the config file.

if impmethod in ["simpute", "npute"]:

    # Call the HAFpipe script with one of the two existing methods.
    rule hafpipe_impute_snp_table:
        input:
            snptable=get_hafpipe_snp_table_dir() + "/{chrom}.csv"
        output:
            # Unnamed output, as this is implicit in HAFpipe Task 2
            get_hafpipe_snp_table_dir() + "/{chrom}.csv" + "." + impmethod
        params:
            tasks="2",
            impmethod=impmethod,
            extra=config["params"]["hafpipe"]["impute-extra"]
        log:
            "logs/hafpipe/impute-" + impmethod + "/{chrom}.log"
        conda:
            "../envs/hafpipe.yaml"
        script:
            "../scripts/hafpipe.py"

elif impmethod != "":

    # Validity check of the custom script.
    if (
        not os.path.exists( config["params"]["hafpipe"]["impute-script"] ) or
        not os.access(config["params"]["hafpipe"]["impute-script"], os.X_OK)
    ):
        raise Exception(
            "User provided impute-script for HAFpipe does not exist or is not executable"
        )

    # Call the user provided script. No need for any of the HAFpipe parameters here.
    rule hafpipe_impute_snp_table:
        input:
            snptable=get_hafpipe_snp_table_dir() + "/{chrom}.csv"
        output:
            # Unnamed output, as this is implicit in the user script
            get_hafpipe_snp_table_dir() + "/{chrom}.csv" + "." + impmethod
        log:
            "logs/hafpipe/impute-" + impmethod + "/{chrom}.log"
        conda:
            # We use the custom conda env if the user provided it,
            # or just re-use the hafpipe env, for simplicity.
            (
                config["params"]["hafpipe"]["impute-conda"]
                if config["params"]["hafpipe"]["impute-conda"] != ""
                else "../envs/hafpipe.yaml"
            )
        shell:
            config["params"]["hafpipe"]["impute-script"] + " {input.snptable}"

# Helper to get the SNP table for a given chromosome. According to the `impmethod` config setting,
# this is either the raw table from Task 1 above, or the imputed table from Task 2, with either one
# of the established methods of HAFpipe, or a custom method/script provided by the user.
def get_hafpipe_snp_table(wildcards):
    base = get_hafpipe_snp_table_dir() + "/" + wildcards.chrom + ".csv"
    if config["params"]["hafpipe"]["impmethod"] in ["", "none"]:
        return base
    else:
        return base + "." + config["params"]["hafpipe"]["impmethod"]

# =================================================================================================
#     HAFpipe Tasks 3 & 4:  Infer haplotype frequencies & Calculate allele frequencies
# =================================================================================================

# We use the merged bam files of all bams per sample.
# The below rule is the same as the mpileup_merge_units rule in pileup.smk, but we do want this
# separate implemenation here, to have a bit more control of where the files go, and to stay
# independent of the mpileup rules. Bit of code duplication, might refactor in the future though.

rule hafpipe_merge_bams:
    input:
        get_sample_bams_wildcards # provided in mapping.smk
    output:
        temp("hafpipe/bam/{sample}.merged.bam")
    params:
        config["params"]["samtools"]["merge"]
    threads:
        config["params"]["samtools"]["merge-threads"]
    log:
        "logs/samtools/hafpipe/merge-{sample}.log"
    wrapper:
        "0.74.0/bio/samtools/merge"

# We run both steps in one rule, just to keep it simple. This requires more runtime for that
# rule, but might be better to save more cluster job submissions? Might refactor to split this.

rule hafpipe_frequencies:
    input:
        bamfile="hafpipe/bam/{sample}.merged.bam",     # provided above
        baifile="hafpipe/bam/{sample}.merged.bam.bai", # provided via bam_index rule in mapping.smk
        snptable=get_hafpipe_snp_table,                # provided above
        refseq=config["data"]["reference"]["genome"]
    output:
        # We currently just specify the output file names here as HAFpipe produces them.
        # Might want to refactor in the future to be able to provide our own names,
        # and have the script rename the files automatically.
        # We don't name the outputs here, as those names are implicit in HAFpipe anyway.
        "hafpipe/frequencies/{sample}.merged.bam.{chrom}.freqs",
        "hafpipe/frequencies/{sample}.merged.bam.{chrom}.afSite"
    params:
        tasks="3,4",
        outdir="hafpipe/frequencies/",
        extra=config["params"]["hafpipe"]["frequencies-extra"]
    log:
        "logs/hafpipe/frequencies/{sample}.{chrom}.log"
    conda:
        "../envs/hafpipe.yaml"
    script:
        "../scripts/hafpipe.py"

# =================================================================================================
#     HAFpipe Collect All
# =================================================================================================

# Get the afSite file list. As this is task 4, task 3 will also be executed,
# but we keep it simple here and only request the final files.
def collect_all_hafpipe_frequencies(wildcards):
    # We use the fai file of the ref genome to cross-check the list of chromosomes in the config.
    # We use a checkpoint to create the fai file from our ref genome, which gives us the chrom names.
    # Snakemake then needs an input function to work with the fai checkpoint here.
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand(
        "hafpipe/frequencies/{sample}.merged.bam.{chrom}.afSite",
        sample=config["global"]["sample-names"],
        chrom=get_hafpipe_chromosomes( fai )
    )

# Simple rule that requests all hafpipe af files, so that they get computed.
# Will probably extend this in the future to a rule that combines all of them into one file.
rule all_hafpipe:
    input:
        collect_all_hafpipe_frequencies

localrules: hafpipe_collect
