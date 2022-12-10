import os
import gzip

# =================================================================================================
#     Setup
# =================================================================================================

# We use a rule for the setup in order to ensure that this is only called once,
# even if there are multiple instances of the downstream roles.

def get_packages_dir():
    # Hard coded paths, within our structure. If we ever want to change the paths where these
    # dependencies are stored, they need to be adjusted in scripts/hafpipe.py as well.
    # Snakemake suggests to use `get_source()`, but it is actually called `source_path()`:
    # https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-does-snakemake-interpret-relative-paths
    # However, neither of them is working... So we use a different way of getting the path
    # to our Snakemake base directory: https://stackoverflow.com/a/73202976/4184258
    return os.path.join( snakemake.workflow.basedir, "packages" )

def get_hafpipe_bins():
    harp_bin    = os.path.join( get_packages_dir(), "harp/bin/harp" )
    hafpipe_bin = os.path.join( get_packages_dir(), "hafpipe/HAFpipe_wrapper.sh" )
    return [ harp_bin, hafpipe_bin ]

rule hafpipe_setup:
    output:
        get_hafpipe_bins()
    params:
        packages_path = get_packages_dir()
    log:
        "logs/hafpipe/setup.log"
    shell:
        "{params.packages_path}/setup-hafpipe.sh >> {log} 2>&1"

localrules:
    hafpipe_setup

# =================================================================================================
#     Input File Check
# =================================================================================================

# Task 1 of HAFpipe (Make SNP Table) expects the founder VCF to contain "PASS" in the FILTER column,
# while the common alternative notation "." does not work. We want to warn users about this,
# as this is just nasty to debug and track down why the SNP table remains empty...
#
# Update: We switched to our our HAF-pipe fork, which fixes that problem already.
# We keep the functions around here for reference, as they might be useful for someone using
# the original HAF-pipe, but we have deactivated their usage below for now.

# Open a file, potentially gzipped, inspired by https://stackoverflow.com/a/16816627
def gzip_opener(filename):
    f = open(filename,'rb')
    if( f.read(2) == b'\x1f\x8b' ):
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f

# Check that the VCF contains PASS rather than . for the FILTER column.
def check_hafpipe_founder_vcf( vcf ):
    with gzip_opener(vcf) as file:
        idx = -1
        dot = 0
        pss = 0
        cnt = 0
        for line in file:
            line = str(line)

            # Skip first header lines
            if line.startswith( "##" ):
                continue
            # Only check first 1k lines. Hope that's enough.
            if cnt > 1000:
                break
            spl = line.split( "\t" )

            # Look at the main header line and find the column of the FILTER
            if line.startswith( "#" ):
                if "FILTER" in spl:
                    idx = spl.index( "FILTER" )
                else:
                    raise Exception("Invalid HAF-pipe founder VCF file without FILTER column")
                continue

            # Look at the genome positions, and count the values found in the FILTER column
            if spl[idx] == ".":
                dot += 1
            elif spl[idx] == "PASS":
                pss += 1
            cnt += 1

        # Check the results
        if dot > pss:
            raise Exception(
                "The FILTER column in the HAF-pipe founder VCF file marks variants with a '.' dot, "
                "typically indicating a missing filter value, and hence likely meaning that the "
                "variant is not filtered out. HAF-pipe however expectes the lines to contain 'PASS' "
                "instead to mark the variants to be used for the SNP table. You hence need to "
                "convert your VCF to use 'PASS' instead of '.' for HAF-pipe to work. "
                "See also https://github.com/petrov-lab/HAFpipe-line/issues/6"
            )

# The below is deactivated, as we switched to our HAF-pipe fork, which fixes that problem already.
# We run this check every time... Bit wasteful, and we could make it a rule that is only
# exectuted once, but that would give the error message in a log file, instead of the main
# Snakemake output. So for now, we keep it like this, so that users get the error immediately.
# if os.path.exists( config["params"]["hafpipe"]["founder-vcf"] ):
#     check_hafpipe_founder_vcf( config["params"]["hafpipe"]["founder-vcf"] )

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
    if len(haf_chrs) == 0:
        haf_chrs = ref_chrs
    haf_chrs = [ str(v) for v in haf_chrs ]
    for chr in haf_chrs:
        if not chr in ref_chrs:
            raise Exception(
                "Chromosome '" + chr + "' specified via the config `params: hafpipe: chromosomes` " +
                "list for running HAFpipe is not part of the reference genome."
            )
    return haf_chrs

# Same function as above, but wrapped to be used in a rule... Snakemake can be complicated...
def get_hafpipe_chromosomes_list(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return get_hafpipe_chromosomes( fai )

# We allow users to specify a directory for the snp table, to avoid recomputation of Task 1.
def get_hafpipe_snp_table_dir():
    cfg_dir = config["params"]["hafpipe"]["snp-table-dir"]
    if cfg_dir:
        return cfg_dir.rstrip('/')
    return "hafpipe/snp-tables"

rule hafpipe_snp_table:
    input:
        vcf=config["params"]["hafpipe"]["founder-vcf"],
        bins=get_hafpipe_bins()
    output:
        snptable   = get_hafpipe_snp_table_dir() + "/{chrom}.csv",
        alleleCts  = get_hafpipe_snp_table_dir() + "/{chrom}.csv.alleleCts",
        # We request both the numeric table as produced by the `numeric_SNPtable.R` script,
        # as well as its compressed binary form. The R script needs an _insane_ amount of memory
        # for larger founder VCF files, but might fail silently when going out of memory.
        # In that case, HAF-pipe just continues as if nothing happend, producing a valid, yet empty,
        # compressed file, which will then lead to errors down the line.
        # So here, by requiring both files, we at least get notified if the R script failed.
        # We further raise an exception in our script when we detect this issue,
        # in order to better inform the user about the situation and how to fix this.
        numeric    = get_hafpipe_snp_table_dir() + "/{chrom}.csv.numeric",
        numericbgz = get_hafpipe_snp_table_dir() + "/{chrom}.csv.numeric.bgz",
        done=touch( get_hafpipe_snp_table_dir() + "/{chrom}.done" )
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

    # Deactivated this warning now, as we switchted to our own HAF-pipe fork, which fixes this.
    # if impmethod == "npute":
        # No comment...
        # logger.warning(
        #     "Using HAF-pipe with SNP table imputation method 'npute' is likely going to fail: "
        #     "We are using Python >= 3.7 in grenepipe, whereas npute requires Pyhon 2.*, "
        #     "and there is unfortunately no easy way to fix this. If you require npute, "
        #     "and get error messages here, please submit an issue to "
        #     "https://github.com/moiexpositoalonsolab/grenepipe/issues, "
        #     "so that we know about this and can try to find a solution.\n"
        #     "Alternatively, you can use the custom impmethod capability of grenepipe "
        #     "to run npute yourself, by providing your own script that runs it.\n"
        # )

    # Call the HAFpipe script with one of the two existing methods.
    rule hafpipe_impute_snp_table:
        input:
            snptable=get_hafpipe_snp_table_dir() + "/{chrom}.csv",
            bins=get_hafpipe_bins()
        output:
            # Unnamed output, as this is implicit in HAFpipe Task 2
            get_hafpipe_snp_table_dir() + "/{chrom}.csv." + impmethod,
            touch(get_hafpipe_snp_table_dir() + "/{chrom}.csv." + impmethod + ".done")
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
            get_hafpipe_snp_table_dir() + "/{chrom}.csv." + impmethod,
            touch(get_hafpipe_snp_table_dir() + "/{chrom}.csv." + impmethod + ".done")
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

# There is an issue with the imputed SNP table being indexed multiple times in parallel environments
# such as clusters, see https://github.com/petrov-lab/HAFpipe-line/issues/5
# We here circumvent this when using imputation, by running the indexing ourselves...
# HAF-pipe is broken, so that's what it takes for now to get this to work properly.
# We use a dummy `flag` file to trigger this step, and ensure that it's executed once per chrom.
# We do not want to require the index files directly in downstream rules (Task 4),
# as they are already created for the base table in the cases without imputation,
# which would confuse snakemake if they already existed...
#
# Update: We fixed this in our HAF-pipe fork, so that this is not longer needed when using
# one of the two HAF-pipe imputation methods. We still need this though for user-defined imputation
# methods.
if impmethod == "":

    # Without imputation, the files should already be there, so we do not need to do anything,
    # just create our indicator trigger file.
    rule hafpipe_snp_table_indices:
        input:
            snptable  = get_hafpipe_snp_table_dir() + "/{chrom}.csv",
            alleleCts = get_hafpipe_snp_table_dir() + "/{chrom}.csv.alleleCts",
            numeric   = get_hafpipe_snp_table_dir() + "/{chrom}.csv.numeric.bgz"
        output:
            flag      = get_hafpipe_snp_table_dir() + "/{chrom}.csv.flag"
        shell:
            "touch {output.flag}"

    localrules:
        hafpipe_snp_table_indices

else:

    # With imputation, we run the scripts that are run by HAF-pipe internally,
    # and mark them as output, to have snakemake validate that they were in fact created.
    rule hafpipe_snp_table_indices:
        input:
            snptable  = get_hafpipe_snp_table_dir() + "/{chrom}.csv." + impmethod,
            bins      = get_hafpipe_bins() # require that HAF-pipe scripts are there
        output:
            alleleCts = get_hafpipe_snp_table_dir() + "/{chrom}.csv." + impmethod + ".alleleCts",
            numeric   = get_hafpipe_snp_table_dir() + "/{chrom}.csv." + impmethod + ".numeric.bgz",
            flag      = get_hafpipe_snp_table_dir() + "/{chrom}.csv." + impmethod + ".flag"
        params:
            hp_scripts = get_packages_dir() + "/hafpipe/scripts"
            # The below code tries to switch where the scripts are found for the old and our
            # new improved version of HAF-pipe... it however does not work if HAF-pipe has not
            # been downloaded yet by the hafpipe_setup rule. We have fixed this now to always
            # use our new version anyway, so we don't need the switch any more.
            # hp_scripts = (
            #     get_packages_dir() + "/hafpipe/scripts"
            #     if os.path.exists( get_packages_dir() + "/hafpipe/scripts" )
            #     else get_packages_dir() + "/hafpipe"
            # )
        log:
            "logs/hafpipe/impute-" + impmethod + "/{chrom}-indices.log"
        conda:
            "../envs/hafpipe.yaml"
        shell:
            # Calls from HAF-pipe, replicated here for our purposes,
            # and with additional checks for whether the files already exist, to not repeat effort.
            # As we are still creating a new file (the flag file), we avoid that snakemake complains
            # if the other two files already exist (because snakemake has a check for this, to
            # avoid that files just "appear" without us intending to have created them).
            "if [ ! -e {input.snptable}{impmethod}.alleleCts ]; then "
            "    echo \"counting alleles in {input.snptable}\" >> {log} 2>&1 ; "
            "    {params.hp_scripts}/count_SNPtable.sh {input.snptable} >> {log} 2>&1 ; "
            "fi ; "
            "if [ ! -e {input.snptable}{impmethod}.numeric.bgz ]; then "
            "    echo \"preparing {input.snptable} for allele frequency calculation\" >> {log} 2>&1 ; "
            "    {params.hp_scripts}/prepare_SNPtable_for_HAFcalc.sh {input.snptable} >> {log} 2>&1 ; "
            "fi ; "
            "touch {output.flag}"

# Helper to get the SNP table for a given chromosome. According to the `impmethod` config setting,
# this is either the raw table from Task 1 above, or the imputed table from Task 2, with either one
# of the established methods of HAFpipe, or a custom method/script provided by the user.
def get_hafpipe_snp_table(wildcards):
    base = get_hafpipe_snp_table_dir() + "/" + wildcards.chrom + ".csv"
    if config["params"]["hafpipe"]["impmethod"] in ["", "none"]:
        return base
    else:
        return base + "." + config["params"]["hafpipe"]["impmethod"]

# We need another helper to process wild cards, requesting the dummy `flag` indicator file
# that ensures that the snp table index files (alleleCt and numeric) are created above.
def get_hafpipe_snp_table_flag(wildcards):
    return get_hafpipe_snp_table(wildcards) + ".flag"

# =================================================================================================
#     HAFpipe Tasks 3 & 4:  Infer haplotype frequencies & Calculate allele frequencies
# =================================================================================================

# We use the merged bam files of all bams per sample.
# The below rule is the same as the mpileup_merge_unit_bams rule in pileup.smk, and in the qualimap
# merging in qc.smk, but we do want this separate implemenation here, to have a bit more control
# of where the files go, and to stay independent of the mpileup rules. Bit of code duplication,
# might refactor in the future though.

rule hafpipe_merge_unit_bams:
    input:
        get_sample_bams_wildcards # provided in mapping.smk
    output:
        # Making the bam files temporary is a bit dangerous, as an error in the calling of Task 4
        # after it has already written the header of the output file will lead to snakemake thinking
        # that the output is valid, hence deleting the bam files... Then, for re-running Task 4 we
        # need to first create the bam files again, which will then update Task 3, which takes ages...
        # But let's assume that all steps work ;-)
        (
            "hafpipe/bam/{sample}.merged.bam"
            if config["params"]["hafpipe"].get("keep-intermediates", True)
            else temp("hafpipe/bam/{sample}.merged.bam")
        ),
        touch("hafpipe/bam/{sample}.merged.done")
    params:
        config["params"]["samtools"]["merge"]
    threads:
        config["params"]["samtools"]["merge-threads"]
    log:
        "logs/samtools/hafpipe/merge-{sample}.log"
    wrapper:
        "0.74.0/bio/samtools/merge"

# We run the two steps separately, so that if Task 4 fails,
# Task 3 does not have to be run again, hence saving time.

rule hafpipe_haplotype_frequencies:
    input:
        bamfile="hafpipe/bam/{sample}.merged.bam",     # provided above
        baifile="hafpipe/bam/{sample}.merged.bam.bai", # provided via bam_index rule in mapping.smk
        snptable=get_hafpipe_snp_table,                # provided above
        flag=get_hafpipe_snp_table_flag,
        refseq=config["data"]["reference-genome"],
        bins=get_hafpipe_bins()
    output:
        # We currently just specify the output file names here as HAFpipe produces them.
        # Might want to refactor in the future to be able to provide our own names,
        # and have the script rename the files automatically.
        freqs=(
            "hafpipe/frequencies/{sample}.merged.bam.{chrom}.freqs"
            if config["params"]["hafpipe"].get("keep-intermediates", True)
            else temp("hafpipe/frequencies/{sample}.merged.bam.{chrom}.freqs")
        ),
        done=touch("hafpipe/frequencies/{sample}.merged.bam.{chrom}.freqs.done")
    params:
        tasks="3",
        outdir="hafpipe/frequencies",
        extra=config["params"]["hafpipe"]["haplotype-frequencies-extra"]
    log:
        "logs/hafpipe/frequencies/haplotype-{sample}.{chrom}.log"
    conda:
        "../envs/hafpipe.yaml"
    script:
        "../scripts/hafpipe.py"

rule hafpipe_allele_frequencies:
    input:
        bamfile="hafpipe/bam/{sample}.merged.bam",     # provided above
        baifile="hafpipe/bam/{sample}.merged.bam.bai", # provided via bam_index rule in mapping.smk
        snptable=get_hafpipe_snp_table,                # provided above
        flag=get_hafpipe_snp_table_flag,
        freqs="hafpipe/frequencies/{sample}.merged.bam.{chrom}.freqs", # from Task 3 above
        bins=get_hafpipe_bins()
    output:
        # Same as above: just expect the file name as produced by HAFpipe.
        afSite=(
            "hafpipe/frequencies/{sample}.merged.bam.{chrom}.afSite"
            if config["params"]["hafpipe"].get("keep-intermediates", True)
            else temp("hafpipe/frequencies/{sample}.merged.bam.{chrom}.afSite")
        ),
        done=touch("hafpipe/frequencies/{sample}.merged.bam.{chrom}.afSite.done")
    params:
        tasks="4",
        outdir="hafpipe/frequencies",
        extra=config["params"]["hafpipe"]["allele-frequencies-extra"]
    log:
        "logs/hafpipe/frequencies/allele-{sample}.{chrom}.log"
    conda:
        "../envs/hafpipe.yaml"
    script:
        "../scripts/hafpipe.py"

# =================================================================================================
#     HAFpipe Concat Sample
# =================================================================================================

# Get the afSite file list for a given sample. As this is the result of task 4,
# task 3 will also be executed, but we keep it simple here and only request the final files.
def collect_sample_hafpipe_allele_frequencies(wildcards):
    # We use the fai file of the ref genome to cross-check the list of chromosomes in the config.
    # We use a checkpoint to create the fai file from our ref genome, which gives us the chrom names.
    # Snakemake then needs an input function to work with the fai checkpoint here.
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand(
        "hafpipe/frequencies/{sample}.merged.bam.{chrom}.afSite",
        sample=wildcards.sample,
        chrom=get_hafpipe_chromosomes( fai )
    )

# Concat all afSite files for a given sample.
# The script assumes the exact naming scheme that we use above, so it is not terribly portable...
rule hafpipe_concat_sample_allele_frequencies:
    input:
        collect_sample_hafpipe_allele_frequencies
    output:
        # This is the file name produced by the script. For now we do not allow to change this.
        table="hafpipe/samples/{sample}.csv" + (
            ".gz" if config["params"]["hafpipe"].get("compress-sample-tables", False) else ""
        ),
        done=touch("hafpipe/samples/{sample}.done")
    params:
        # The rule needs access to the list of chromosomes, and to the sample.
        sample="{sample}",
        chroms=get_hafpipe_chromosomes_list,

        # We might want to compress the output per sample.
        compress=config["params"]["hafpipe"].get("compress-sample-tables", False)
    log:
        "logs/hafpipe/samples/{sample}.log"
    script:
        "../scripts/hafpipe-concat.py"

# Simply request all the above sample files.
rule hafpipe_collect_concat_samples:
    input:
        tables=expand(
            "hafpipe/samples/{sample}.csv" + (
                ".gz" if config["params"]["hafpipe"].get("compress-sample-tables", False) else ""
            ),
            sample=config["global"]["sample-names"]
        )
    output:
        done=touch("hafpipe/samples.done")

localrules: hafpipe_collect_concat_samples

# =================================================================================================
#     HAFpipe Merge All
# =================================================================================================

# Get the afSite file list. As this is the result of task 4, task 3 will also be executed,
# but we keep it simple here and only request the final files.
def collect_all_hafpipe_allele_frequencies(wildcards):
    # We use the fai file of the ref genome to cross-check the list of chromosomes in the config.
    # We use a checkpoint to create the fai file from our ref genome, which gives us the chrom names.
    # Snakemake then needs an input function to work with the fai checkpoint here.
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand(
        "hafpipe/frequencies/{sample}.merged.bam.{chrom}.afSite",
        sample=config["global"]["sample-names"],
        chrom=get_hafpipe_chromosomes( fai )
    )

# Merge all afSite files produced above, for all samples and all chromsomes.
# The script assumes the exact naming scheme that we use above, so it is not terribly portable...
rule hafpipe_merge_allele_frequencies:
    input:
        # We only request the input files here so that snakemake knows that we need them.
        # The script that we are running does not access the input list at all, as it is just a
        # big mess of files. We instead want the structured way of finding our files using the
        # lists of samples and chromosomes that we hand over via the params below.
        # This is unfortunately necessary, as HAFpipe afSite files do not contain any information
        # on their origins (samples and chromosomes) other than their file names, so we have to
        # work with that... See the script for details.
        collect_all_hafpipe_allele_frequencies
    output:
        # This is the file name produced by the script. For now we do not allow to change this.
        table="hafpipe/all.csv" + (
            ".gz" if config["params"]["hafpipe"].get("compress-merged-table", False) else ""
        ),
        done=touch("hafpipe/all.done")
    params:
        # We are potentially dealing with tons of files, and cannot open all of them at the same
        # time, due to OS limitations, check `ulimit -n` for example. When this param is set to 0,
        # we try to use that upper limit (minus some tolerance). However, if that fails, this value
        # can be manually set to a value below the limit, to make it work, e.g., 500
        concurrent_files=0,

        # The rule needs access to lists of samples and chromosomes,
        # and we give it the base path of the afSite files as well, for a little bit of flexibility.
        samples=config["global"]["sample-names"],
        chroms=get_hafpipe_chromosomes_list,

        # We provide the paths to the input directory here. The output will be written to the
        # parent directory of that (so, to "hafpipe").
        # Ugly, but we are dealing with HAFpipe uglines here, and that seems to be easiest for now.
        base_path="hafpipe/frequencies",

        # We might want to compress the final output.
        compress=config["params"]["hafpipe"].get("compress-merged-table", False)
    log:
        "logs/hafpipe/merge-all.log"
    script:
        "../scripts/hafpipe-merge.py"

# =================================================================================================
#     HAFpipe All
# =================================================================================================

# Simply request all the above afSite files.
# We always request to run this, so that HAF-pipe is run even when the two make table config
# options are not set. This is so that users who just want HAF-pipe out as-is are able to get that.
rule hafpipe_collect_allele_frequencies:
    input:
        collect_all_hafpipe_allele_frequencies
    output:
        done=touch("hafpipe/afSite.done")

localrules: hafpipe_collect_allele_frequencies

# Simple rule that requests all hafpipe af files, so that they get computed.
# This requests the completely merged table, and/or the per-sample concatenated tables,
# and always the basic HAF-pipe output, so that those are always produced independently of the tables.
rule all_hafpipe:
    input:
        "hafpipe/afSite.done",
        "hafpipe/all.csv" + (
            ".gz" if config["params"]["hafpipe"].get("compress-merged-table", False) else ""
        ) if config["params"]["hafpipe"].get("make-merged-table", False) else [],
        "hafpipe/samples.done" if config["params"]["hafpipe"].get("make-sample-tables", False) else []

localrules: all_hafpipe
