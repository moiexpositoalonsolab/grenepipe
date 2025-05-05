# =================================================================================================
#     Make Contig Shards
# =================================================================================================

# Special rule for the case that we want to split the contigs for the genomics db import,
# but are not using contig groups. In that case, we could pretend to have contig groups,
# and create fake ones that each contain a single chromosome or contig of the reference.
# However, that would mean remodelling some other code here as well which instead expects
# proper contig names, such as the `call_variants` rule based on get_gatk_intervals().
# So insted of refactoring those, we instead here create a bed file for the chromosomes,
# without pretending this to be contig groups.
rule make_contig_bed:
    input:
        fai=get_fai,
    output:
        bed="calling/contig-shards/{contig}/contig.bed"
    params:
        contig="{contig}"
    run:
        # Produce only the matching contig and write 0–length interval
        contig_lengths = get_contig_lengths(input.fai)
        if params.contig not in contig_lengths:
            raise Exception("Invalid contig not in fai file: " + str(params.contig))
        with open(output.bed, "w") as out:
            name = params.contig
            length = contig_lengths[params.contig]
            out.write(f"{name}\t0\t{length}\n")


localrules:
    make_contig_bed,


# Checkpoint: tile each contig into fixed-width intervals (with optional padding).
# For some reason, GATK _needs_ the --interval-merging-rule option to be set to OVERLAPPING_ONLY.
# Apparently, not setting this leads to an error asking to set it.
# Why that then is an option at all is beyond my comprehension. GATK, WTF.
checkpoint preprocess_contig_shard:
    input:
        dict = genome_dict(),
        ref  = config["data"]["reference-genome"],
        # contigs=contigs_groups_input,
        contig = (
            "calling/contig-groups/{contig}.bed"
            if config["settings"].get("contig-group-size", 0) > 0
            else "calling/contig-shards/{contig}/contig.bed"
        )
    output:
        interval_list = "calling/contig-shards/{contig}/contig.interval_list",
        shard_list = "calling/contig-shards/{contig}/shards.interval_list",
    params:
        bin_length = config["params"]["gatk"]["GenomicsDBImport-interval-size"],
        padding    = config["params"]["gatk"]["GenomicsDBImport-interval-padding"],
    log:
        "logs/calling/preprocess-contig-shard/{contig}.log",
    conda:
        "../envs/gatk.yaml"
    shell:
        """
        # convert to interval_list
        mkdir -p $(dirname {output.interval_list})
        gatk BedToIntervalList \
            -I {input.contig} \
            -SD {input.dict} \
            -O {output.interval_list} \
            &>> {log}
        # tile into fixed-size bins
        mkdir -p $(dirname {output.shard_list})
        gatk PreprocessIntervals \
            -R {input.ref} \
            -L {output.interval_list} \
            --bin-length {params.bin_length} \
            --padding {params.padding} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O {output.shard_list} \
            &>> {log}
        """


localrules:
    preprocess_contig_shard,


# Extract one interval per file
rule extract_contig_shard:
    input:
        # shard_list=checkpoints.preprocess_contig_shard.output.shard_list
        shard_list="calling/contig-shards/{contig}/shards.interval_list"
    output:
        shard = "calling/contig-shards/{contig}/shard-{shard}.interval_list"
    run:
        header, data = [], []
        for line in open(input.shard_list):
            if line.startswith('@'):
                header.append(line)
            else:
                data.append(line)
        i = int(wildcards.shard)
        with open(output.shard, 'w') as out:
            out.writelines(header + [data[i]])


localrules:
    extract_contig_shard,


# helper to list shards for a contig
def get_shard_indices(wc):
    r = checkpoints.preprocess_contig_shard.get(contig=wc.contig)
    # skip header lines
    shards = [l for l in open(r.output.shard_list) if not l.startswith('@')]
    return [i for i in range(len(shards))]
    # return [{"shard": i} for i in range(len(shards))]


# =================================================================================================
#     Genomics DB
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
        # Here we specify the intervals to process, which are defined by a shard.
        # This is an input file from above, making sure that this is created beforehand.
        intervals="calling/contig-shards/{contig}/shard-{shard}.interval_list",
    output:
        db=directory("calling/genomics_db/{contig}/db-shard-{shard}"),
        done=touch("calling/genomics_db/{contig}/db-shard-{shard}.done"),
    params:
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
        "logs/calling/gatk-genomicsdbimport/{contig}/db-shard-{shard}.log",
    benchmark:
        "benchmarks/calling/gatk-genomicsdbimport/{contig}/db-shard-{shard}.log"
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
        genomicsdb="calling/genomics_db/{contig}/db-shard-{shard}",
        genomicsdb_done="calling/genomics_db/{contig}/db-shard-{shard}.done",
        # db = rules.genomics_db_import_shard.output.db,
        # If known variants are set in the config, use them, and require the index file as well.
        known=config["data"]["known-variants"],
        knownidx=(
            config["data"]["known-variants"] + ".tbi" if config["data"]["known-variants"] else []
        ),
        # Same as above, we specify the intervals to process, which are defined by a shard.
        # This is an input file from above, making sure that this is created beforehand.
        intervals="calling/contig-shards/{contig}/shard-{shard}.interval_list",
    output:
        vcf=(
            "calling/genotyped/{contig}/shard-{shard}.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/genotyped/{contig}/shard-{shard}.vcf.gz")
        ),
        done=touch("calling/genotyped/{contig}/shard-{shard}.vcf.gz.done"),
    params:
        extra=" --sequence-dictionary "
        + genome_dict()
        + " "
        + config["params"]["gatk"]["GenotypeGVCFs-extra"],
        java_opts=config["params"]["gatk"]["GenotypeGVCFs-java-opts"],
    threads:
        get_rule_threads("genotype_variants")
    log:
        "logs/calling/gatk-genotype-gvcfs/{contig}/shard-{shard}.log",
    benchmark:
        "benchmarks/calling/gatk-genotype-gvcfs/{contig}/shard-{shard}.log"
    # group:
    #     "gatk_calls_combine"
    conda:
        "../envs/gatk.yaml"
    wrapper:
        "v5.7.0/bio/gatk/genotypegvcfs"


# =================================================================================================
#     Gather Shards
# =================================================================================================


# def merge_vcfs_vcfs_input(wildcards):
#     # r = checkpoints.preprocess_contig_shard.get(contig=wildcards.contig)
#     return expand(
#         "calling/genotyped/{{contig}}/shard-{shard}.vcf.gz",
#         shard=get_shard_indices(wildcards)
#     )
#
# def merge_vcfs_done_input(wildcards):
#     # r = checkpoints.preprocess_contig_shard.get(contig=wildcards.contig)
#     return expand(
#         "calling/genotyped/{{contig}}/shard-{shard}.vcf.gz.done",
#         shard=get_shard_indices(wildcards)
#     )

def merge_vcfs_vcfs_input(wc):
    cp = checkpoints.preprocess_contig_shard.get(**wc)
    # cp = checkpoints.preprocess_contig_shard.get(contig=wc.contig)
    # now that the checkpoint has run, cp.output.shard_list exists
    with open(cp.output.shard_list) as f:
        shards = [l for l in f if not l.startswith("@")]
    return expand(
        "calling/genotyped/{contig}/shard-{shard}.vcf.gz",
        contig=wc.contig,
        shard=list(range(len(shards)))
    )

def merge_vcfs_done_input(wc):
    cp = checkpoints.preprocess_contig_shard.get(**wc)
    # cp = checkpoints.preprocess_contig_shard.get(contig=wc.contig)
    with open(cp.output.shard_list) as f:
        shards = [l for l in f if not l.startswith("@")]
    return expand(
        "calling/genotyped/{contig}/shard-{shard}.vcf.gz.done",
        contig=wc.contig,
        shard=list(range(len(shards)))
    )

rule merge_shard_vcfs:
    input:
        dict = genome_dict(),
        ref  = config["data"]["reference-genome"],
        contigs=contigs_groups_input,
        # vcfs = expand("calling/genotyped/{contig}/shard-{shard}.vcf.gz", get_shard_indices),
        # done = expand("calling/genotyped/{contig}/shard-{shard}.vcf.gz.done", get_shard_indices),
        vcfs=merge_vcfs_vcfs_input,
        done=merge_vcfs_done_input,
    output:
        vcf=(
            "calling/genotyped/all.{contig}.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/genotyped/all.{contig}.vcf.gz")
        ),
        done=touch("calling/genotyped/all.{contig}.vcf.gz.done"),
    params:
        # See duplicates-picard.smk for the reason whe need this on MacOS.
        java_opts=config["params"]["picard"].get("MergeVcfs-java-opts", ""),
        extra=(
            " --USE_JDK_DEFLATER true --USE_JDK_INFLATER true"
            if platform.system() == "Darwin"
            else ""
        ),
    threads:
        get_rule_threads("merge_shard_vcfs")
    log:
        "logs/calling/picard-merge-vcfs/{contig}.log",
    benchmark:
        "benchmarks/calling/picard-merge-vcfs/{contig}.log"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "v5.7.0/bio/picard/mergevcfs"
