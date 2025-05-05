# =================================================================================================
#     Combining Calls
# =================================================================================================


# Recommended way of GATK to combine GVCFs these days.
rule genomics_db_import:
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
        # Same as in the haplotype caller, we need a dummy for the intervals
        # to ensure the files are present.
        intervals_dummy=get_gatk_interval_files,
    output:
        db=directory("calling/genomics_db/{contig}"),
        done=touch("calling/genomics_db/{contig}.done"),
    params:
        # Here, we actually use the intervals to provide them to the wrapper.
        intervals=get_gatk_intervals,
        db_action="create",
        extra=" --reference "
        + config["data"]["reference-genome"]
        + " --sequence-dictionary "
        + genome_dict()
        + " "
        + config["params"]["gatk"].get("GenomicsDBImport-extra", ""),
        java_opts=config["params"]["gatk"].get("GenomicsDBImport-java-opts", ""),
    threads: get_rule_threads("genomics_db_import")
    log:
        "logs/calling/gatk-genomicsdbimport/{contig}.log",
    benchmark:
        "benchmarks/calling/gatk-genomicsdbimport/{contig}.log"
    resources:
        tmpdir=config["params"]["gatk"].get("GenomicsDBImport-temp-dir", ""),
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/genomicsdbimport"


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
        # Get the GenomicsDB input
        genomicsdb="calling/genomics_db/{contig}",
        genomicsdb_done="calling/genomics_db/{contig}.done",
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
