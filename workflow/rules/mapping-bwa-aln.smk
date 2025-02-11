import platform

# =================================================================================================
#     Read Mapping
# =================================================================================================


# We need a special rule that gives us exactly one fastq file. Our get_trimmed_reads() returns
# a list of either one or two fastq files, depending on whether it is a single or paired end
# sample. bwa aln cannot work with this, so we need to call it independently, and hence need
# to split this list accordingly.
def get_single_trimmed_read(wildcards):
    # Get the list that comes from the trimming tool.

    list = get_trimmed_reads(wildcards)
    if len(list) not in [1, 2]:
        raise Exception("Invalid trimming tool result that is neither se nor pe")

    # Return one value of the two, with a bit of internal error checking.
    # Probably more error checking than necessary, but better safe than sorry.
    if wildcards.pair == "1":
        return list[0]
    elif wildcards.pair == "2":
        if len(list) != 2:
            raise Exception("Invalid pair number 2 requested for single end sample")
        return list[1]
    else:
        raise Exception("Invalid pair number that is not in [1,2] requested")


rule map_reads:
    input:
        reads=get_single_trimmed_read,
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        reads_done=get_trimmed_reads_done,
    output:
        (
            "mapping/sai/{sample}-{unit}-{pair}.sai"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapping/sai/{sample}-{unit}-{pair}.sai")
        ),
        touch("mapping/sai/{sample}-{unit}-{pair}.sai.done"),
        # Absolutely no idea why the following does not work. Snakemake complains that the pipe
        # is not used at all, which is simply wrong, as it is clearly used by bwa_sai_to_bam,
        # and also why would this rule here be called if its output file was not requested at all?!
        # Snakemake... no idea what you are doing there again.
        # pipe( "mapping/sai/{sample}-{unit}-{pair}.sai" )
    params:
        index=config["data"]["reference-genome"],
        extra=config["params"]["bwaaln"]["extra"],
    group:
        "mapping"
    log:
        "logs/mapping/bwa-aln/{sample}-{unit}-{pair}.log",
    benchmark:
        "benchmarks/mapping/bwa-aln/{sample}-{unit}-{pair}.log"
    threads: config["params"]["bwaaln"]["threads"]
    # resources:
    # Increase time limit in factors of 2h, if the job fails due to time limit.
    # time = lambda wildcards, input, threads, attempt: int(120 * int(attempt))
    conda:
        "../envs/bwa.yaml"
    wrapper:
        "0.74.0/bio/bwa/aln"


# =================================================================================================
#     Convert to bam
# =================================================================================================

# Apprently, samtools does not create the tmp dir correclty, so we need to take care of this...
if len(config["params"]["samtools"]["temp-dir"]) > 0:
    os.makedirs(config["params"]["samtools"]["temp-dir"], exist_ok=True)

# We use the BWA SAM(SE/PE) wrapper, which can handle both se and pe files at the same time.
# This is really convenient, as this means we don't have to deal with this ourselves. Nice!
# Apparently, it uses the input fastq files to determine if we have se or pe data. Smart.
# Still, we ned a bit of trickery to get it to work.


# Adapted from get_trimmed_reads, but replace fastq by our sai files produced in the rule above.
def get_sai(wildcards):
    if is_single_end(**wildcards):
        # Single end sample.
        return ["mapping/sai/{sample}-{unit}-1.sai".format(**wildcards)]
    elif config["settings"]["merge-paired-end-reads"]:
        # Merged paired-end samples.
        # Here, we rely on the fact that get_trimmed_reads() only returns a single sample for
        # merged paired-end samples, which is then the only one in the list returned from that
        # function. So then, we pretend that this is the "pair == 1" sample, so that the bwa aln
        # rule maps that file.
        return ["mapping/sai/{sample}-{unit}-1.sai".format(**wildcards)]
    else:
        # Paired-end sample.
        return expand(
            "mapping/sai/{sample}-{unit}-{pair}.sai",
            sample=wildcards.sample,
            unit=wildcards.unit,
            pair=[1, 2],
        )

def get_sai_done(wildcards):
    return get_sai(wildcards) + ".done"

def get_bwa_aln_extra(wildcards):
    # We need the read group tags, including `ID` and `SM`, as downstream tools use these.
    # Contrary to bwa mem, this is here specified with lowercase -r, instead of uppercase -R.
    # As if bioinformatics tools were ever consistent...
    rg_tags = "\\t".join(get_read_group_tags(wildcards))
    extra = "-r '@RG\\t" + rg_tags + "' " + config["params"]["bwaaln"]["extra-sam"]
    return extra


rule bwa_sai_to_bam:
    input:
        fastq=get_trimmed_reads,
        sai=get_sai,
        done=get_sai_done,
        ref=config["data"]["reference-genome"],
        # Somehow, the wrapper expects the index extensions to be given,
        # instead of the underlying fasta file... Well, so let's do that.
        # We provide the fasta above as well; it's not used,
        # but might be important as a rule dependency so that it is present.
        idx=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
    output:
        pipe("mapping/sorted/{sample}-{unit}-unclean.bam"),
        # touch("mapping/sorted/{sample}-{unit}-unclean.bam.done"),
    params:
        extra=get_bwa_aln_extra,
        # Sort as we need it.
        sort="samtools",
        sort_order="coordinate",
        sort_extra=config["params"]["samtools"]["sort"],
        tmp_dir=config["params"]["samtools"]["temp-dir"],
    group:
        "mapping"
    log:
        "logs/mapping/bwa-sam/{sample}-{unit}.log",
    benchmark:
        "benchmarks/mapping/bwa-sam/{sample}-{unit}.log"
    conda:
        # The wrapper does not include numpy and pandas as dependencies, but somehow needs them...
        # So we just re-use our normal bwa env, which also workes for the above rule.
        "../envs/bwa.yaml"
    script:
        # We use our own version of the wrapper here, as that wrapper misses temp dirs for
        # samtools sort, causing all kinds of trouble...
        # wrapper:
        #     "0.74.0/bio/bwa/samxe"
        "../scripts/bwa-samxe.py"


# Apparently, yet another bioinformatics tool fail is at play here. The bam files written above
# lead to a SAM validation error: ERROR::INVALID_MAPPING_QUALITY, MAPQ should be 0 for unmapped read
# when opened with Picard MarkDuplicates, see also https://www.biostars.org/p/55830/
# We hence here use Picard CleanSam to clean them up again, so that downstream Picard MarkDuplicates
# can work with those files again.
rule bwa_bam_clean:
    input:
        "mapping/sorted/{sample}-{unit}-unclean.bam",
    output:
        (
            "mapping/sorted/{sample}-{unit}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("mapping/sorted/{sample}-{unit}.bam")
        ),
        touch("mapping/sorted/{sample}-{unit}.bam.done"),
    params:
        # See duplicates-picard.smk for the reason whe need this on MacOS.
        extra=(
            " --USE_JDK_DEFLATER true --USE_JDK_INFLATER true" if platform.system() == "Darwin" else ""
        ),
    group:
        "mapping"
    log:
        "logs/mapping/picard-cleansam/{sample}-{unit}.log",
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CleanSam --INPUT {input[0]} --OUTPUT {output[0]} {params.extra} &> {log}"
        # Somehow, Picard has several active versions with different command line interfaces.
        # Lets' hope that we picked the one that works...
        # "picard CleanSam --INPUT {input} --OUTPUT {output}"
