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
    if len(list) not in [ 1, 2 ]:
        raise Exception( "Invalid trimming tool result that is neither se nor pe" )

    # Return one value of the two, with a bit of internal error checking.
    # Probably more error checking than necessary, but better safe than sorry.
    if wildcards.pair == "1":
        return list[0]
    elif wildcards.pair == "2":
        if len(list) != 2:
            raise Exception( "Invalid pair number 2 requested for single end sample" )
        return list[1]
    else:
        raise Exception( "Invalid pair number that is not in [1,2] requested" )

rule map_reads:
    input:
        reads=get_single_trimmed_read,
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        )
    output:
        "sai/{sample}-{unit}-{pair}.sai"
        # Absolutely no idea why the following does not work. Snakemake complains that the pipe
        # is not used at all, which is simply wrong, as it is clearly used by bwa_sai_to_bam,
        # and also why would this rul here be called if its output file was not requested at all?!
        # Snakemake... no idea what you are doing there again.
        # pipe( "sai/{sample}-{unit}-{pair}.sai" )
    params:
        index=config["data"]["reference"]["genome"],
        extra=config["params"]["bwaaln"]["extra"]
    group:
        "mapping"
    log:
        "logs/bwa-aln/{sample}-{unit}-{pair}.log"
    benchmark:
        "benchmarks/bwa-aln/{sample}-{unit}-{pair}.bench.log"
    threads:
        config["params"]["bwaaln"]["threads"]
    resources:
        # Increase time limit in factors of 2h, if the job fails due to time limit.
        time = lambda wildcards, input, threads, attempt: int(120 * int(attempt))
    shadow:
        # See the bwa mem rule for the rationale to use a full shadow directory here.
        "full"
    conda:
        "../envs/bwa.yaml"
    wrapper:
        "0.74.0/bio/bwa/aln"

# =================================================================================================
#     Convert to bam
# =================================================================================================

# We use the BWA SAM(SE/PE) wrapper, which can handle both se and pe files at the same time.
# This is really convenient, as this means we don't have to deal with this ourselves. Nice!
# Apparently, it uses the input fastq files to determine if we have se or pe data. Smart.
# Still, we ned a bit of trickery to get it to work.

# Adapted from get_trimmed_reads, but replace fastq by our sai files produced in the rule above.
def get_sai(wildcards):
    if is_single_end(**wildcards):
        # Single end sample.
        return [ "sai/{sample}-{unit}-1.sai".format(**wildcards) ]
    elif config["settings"]["merge-paired-end-reads"]:
        # Merged paired-end samples.
        # Here, we rely on the fact that get_trimmed_reads() only returns a single sample for
        # merged paired-end samples, which is then the only one in the list returned from that
        # function. So then, we pretend that this is the "pair == 1" sample, so that the bwa aln
        # rule maps that file.
        return [ "sai/{sample}-{unit}-1.sai".format(**wildcards) ]
    else:
        # Paired-end sample.
        return expand(
            "sai/{sample}-{unit}-{pair}.sai",
            sample=wildcards.sample, unit=wildcards.unit, pair=[1, 2]
        )

rule bwa_sai_to_bam:
    input:
        fastq=get_trimmed_reads,
        sai=get_sai,
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        )
    output:
        pipe( "mapped/{sample}-{unit}.sorted-unclean.bam" )
    params:
        index=config["data"]["reference"]["genome"],

        # We need the read group tags, including `ID` and `SM`, as downstream tools use these.
        # Contrary to bwa mem, this is here specified with lowercase -r, instead of uppercase -R.
        # As if bioinformatics tools were ever consistent...
        extra=r"-r '@RG\tID:{sample}\tSM:{sample}' " + config["params"]["bwaaln"]["extra-sam"],

        # Sort as we need it.
        sort="samtools",
        sort_order="coordinate",
        sort_extra=config["params"]["bwaaln"]["extra-sort"]
    group:
        "mapping"
    log:
        "logs/bwa-sam/{sample}-{unit}.log"
    benchmark:
        "benchmarks/bwa-sam/{sample}-{unit}.bench.log"
    conda:
        # The wrapper does not include numpy and pandas as dependencies, but somehow needs them...
        # So we just re-use our normal bwa env, which also workes for the above rule.
        "../envs/bwa.yaml"
    wrapper:
        "0.74.0/bio/bwa/samxe"

# Apparently, yet another bioinformatics tool fail is at play here. The bam files written above
# lead to a SAM validation error: ERROR::INVALID_MAPPING_QUALITY, MAPQ should be 0 for unmapped read
# when opened with Picard MarkDuplicates, see also https://www.biostars.org/p/55830/
# We hence here use Picard CleanSam to clean them up again, so that downstream Picard MarkDuplicates
# can work with those files again.
rule bwa_bam_clean:
    input:
        "mapped/{sample}-{unit}.sorted-unclean.bam"
    output:
        "mapped/{sample}-{unit}.sorted.bam"
    group:
        "mapping"
    log:
        "logs/picard/cleansam/{sample}-{unit}.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard CleanSam INPUT={input} OUTPUT={output} &> {log}"
        # Somehow, Picard has several active versions with different command line interfaces.
        # Lets' hope that we picked the one that works...
        # "picard CleanSam --INPUT {input} --OUTPUT {output}"
