# =================================================================================================
#     Combining Calls
# =================================================================================================


# Old way, using GATK CombineGVCFs, which is slow when run on many samples.
# We still offer it for compatibility and completeness, but recommend using GenomicsDBImport instead.
rule combine_calls:
    input:
        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        refdict=genome_dict(),
        # Get the sample data, including indices.
        gvcfs=expand(
            "calling/called/{sample}-{{contig}}.g.vcf.gz", sample=config["global"]["sample-names"]
        ),
        indices=expand(
            "calling/called/{sample}-{{contig}}.g.vcf.gz.tbi",
            sample=config["global"]["sample-names"],
        ),
        done=expand(
            "calling/called/{sample}-{{contig}}.g.vcf.gz.done",
            sample=config["global"]["sample-names"],
        ),
    output:
        gvcf=(
            "calling/combined/all.{contig}.g.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/combined/all.{contig}.g.vcf.gz")
        ),
        done=touch("calling/combined/all.{contig}.g.vcf.gz.done"),
    params:
        extra=config["params"]["gatk"]["CombineGVCFs-extra"]
        + (
            " --dbsnp " + config["data"]["known-variants"] + " "
            if config["data"]["known-variants"]
            else ""
        ),
        java_opts=config["params"]["gatk"]["CombineGVCFs-java-opts"],
    threads: get_rule_threads("combine_calls")
    log:
        "logs/calling/gatk-combine-gvcfs/{contig}.log",
    benchmark:
        "benchmarks/calling/gatk-combine-gvcfs/{contig}.log"
    # group:
    #     "gatk_calls_combine"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/combinegvcfs"


# =================================================================================================
#     Genotype Variants
# =================================================================================================


rule genotype_variants:
    input:
        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        refdict=genome_dict(),
        # Get the GVCF or GenomicsDB input, depending on which tool is requested in the config.
        gvcf="calling/combined/all.{contig}.g.vcf.gz",
        gvcf_done="calling/combined/all.{contig}.g.vcf.gz.done",
        # If known variants are set in the config, use them, and require the index file as well.
        known=config["data"]["known-variants"],
        knownidx=(
            config["data"]["known-variants"] + ".tbi" if config["data"]["known-variants"] else []
        ),
        # Same as above, we need a dummy for the intervals to ensure the files are present.
        intervals_dummy=get_gatk_interval_files,
    output:
        vcf=(
            "calling/genotyped/all.{contig}.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/genotyped/all.{contig}.vcf.gz")
        ),
        done=touch("calling/genotyped/all.{contig}.vcf.gz.done"),
    params:
        # Again, we here use the intervals to provide them to the wrapper.
        intervals=get_gatk_intervals,
        extra=" --sequence-dictionary "
        + genome_dict()
        + " "
        + config["params"]["gatk"]["GenotypeGVCFs-extra"],
        java_opts=config["params"]["gatk"]["GenotypeGVCFs-java-opts"],
    threads: get_rule_threads("genotype_variants")
    log:
        "logs/calling/gatk-genotype-gvcfs/{contig}.log",
    benchmark:
        "benchmarks/calling/gatk-genotype-gvcfs/{contig}.log"
    # group:
    #     "gatk_calls_combine"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/genotypegvcfs"
