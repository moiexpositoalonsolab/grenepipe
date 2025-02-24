# =================================================================================================
#     Variant Calling
# =================================================================================================


rule call_variants:
    input:
        # Need the ref genome, as well as its indices.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        # Get the sample data.
        samples=get_sample_bams_wildcards,
        indices=get_sample_bais_wildcards,
        done=get_sample_bams_wildcards_done,
        # If we use restricted regions, set them here. If not, empty, which will propagate to the
        # get_mpileup_params function as well. Same for small contig groups.
        # regions="calling/regions/{contig}.bed" if config["settings"].get("restrict-regions") else []
        regions=(
            "calling/regions/{contig}.bed"
            if (config["settings"].get("restrict-regions"))
            else (
                "calling/contig-groups/{contig}.bed"
                if (config["settings"].get("contig-group-size"))
                else []
            )
        ),
    output:
        gvcf=(
            "calling/called/{sample}.{contig}.g.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/called/{sample}.{contig}.g.vcf.gz")
        ),
        gtbi="calling/called/{sample}.{contig}.g.vcf.gz.tbi",
        done=touch("calling/called/{sample}.{contig}.g.vcf.gz.done"),
    params:
        # Optional parameters for bcftools mpileup (except -g, -f).
        mpileup=get_mpileup_params,
        # Optional parameters for bcftools call (except -v, -o, -m).
        call=config["params"]["bcftools"]["call"],
    log:
        "logs/calling/bcftools-call/{sample}.{contig}.log",
    benchmark:
        "benchmarks/calling/bcftools-call/{sample}.{contig}.log"
    conda:
        "../envs/bcftools.yaml"
    threads: config["params"]["bcftools"]["threads"]
    shell:
        # We need a normalization step in between, because somehow bcftools spits out data
        # in some non-normal format by default... see https://www.biostars.org/p/404061/#404245
        # See https://www.biostars.org/p/397616/ for details on the annotations
        # "    -a \"AD,ADF,ADR,DP,SP,INFO/AD,INFO/ADF,INFO/ADR,FORMAT/AD,FORMAT/DP\" --gvcf 0 | "
        "("
        "bcftools mpileup "
        "    {params.mpileup} --fasta-ref {input.ref} --output-type u {input.samples} "
        "    -a 'INFO/AD' -a 'FORMAT/AD' -a 'FORMAT/DP' --gvcf 0 | "
        "bcftools call "
        "    --threads {threads} {params.call} --gvcf 0 --output-type u | "
        "bcftools norm "
        "    --fasta-ref {input.ref} "
        "    --output-type z -o {output.gvcf} ; "
        "bcftools index --tbi {output.gvcf}"
        ") &> {log}"


# Potential reference implementation:
# https://gist.github.com/lindenb/f5311887f5a5df0abfaf487df4ae2858

# =================================================================================================
#     Combining Contigs
# =================================================================================================


rule combine_contig:
    input:
        # Need the ref genome, as well as its indices.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
        # Get the sample data, including indices, which are produced above.
        gvcfs=expand(
            "calling/called/{sample}.{{contig}}.g.vcf.gz", sample=config["global"]["sample-names"]
        ),
        indices=expand(
            "calling/called/{sample}.{{contig}}.g.vcf.gz.tbi",
            sample=config["global"]["sample-names"],
        ),
        done=expand(
            "calling/called/{sample}.{{contig}}.g.vcf.gz.done",
            sample=config["global"]["sample-names"],
        ),
    output:
        gvcf="calling/combined/all.{contig}.g.vcf.gz",
        gtbi="calling/combined/all.{contig}.g.vcf.gz.tbi",
        gvcflist="calling/combined/all.{contig}.g.txt",
        done=touch("calling/combined/all.{contig}.g.vcf.gz.done"),
    log:
        "logs/calling/bcftools-combine-contig/{contig}.log",
    benchmark:
        "benchmarks/calling/bcftools-combine-contig/{contig}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        # Store the input gvcfs in a list file to be accessed by bcftool; we are using the same
        # order as the samples table, to keep the gvcf headers identical (needed for concat later).
        # Let's hope that this works with large numbers of files without exceeding the ulimit
        # for number of files that can be kept open. Needs to be tested on a large dataset though...
        "echo {input.gvcfs} | sed 's/ /\\n/g' > {output.gvcflist} ; "
        "("
        "bcftools merge --merge all --gvcf {input.ref} --file-list {output.gvcflist} --output-type u | "
        "bcftools convert --gvcf2vcf --fasta-ref {input.ref} --output-type u | "
        "bcftools view -m 2 --output-type u | "
        "bcftools norm --check-ref ws --fasta-ref {input.ref} --output-type z -o {output.gvcf} ; "
        "bcftools index --tbi {output.gvcf}"
        ") &> {log}"


# =================================================================================================
#     Combining All
# =================================================================================================


# We use the fai file of the reference genome to obtain a list of contigs. This is used as input in
# the rule below, which then request all per-contig gvcfs. We need to provide this as an input
# function for the rule input, so that the fai snakemake checkpoint gets executed before evaluating.
# The list of contigs is either just the list of chromosomes, or, if we use contig grouping in the
# config file, this is the list of the group names.
def combined_contig_gvcfs(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("calling/combined/all.{contig}.g.vcf.gz", contig=get_contigs(fai))


# Also need the done files to make sure snakemake doesn't mess this up.
def combined_contig_done(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("calling/combined/all.{contig}.g.vcf.gz.done", contig=get_contigs(fai))


# We also need a comma-separated list of the contigs, so that bcftools can output
# the concatenated entries in the correct order as given in the fai.
# For this, we use the same technique of using the fai checkpoint as before.
# Here however, we use the actual list of chromosome names, and never the list of contig group names.
def combined_contig_order(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    contig_list = []
    with open(fai, "r") as f:
        for line in f:
            contig = line.split("\t")[0].strip()
            contig_list.append(contig)
    return ",".join(contig_list)


rule combine_all:
    input:
        # fai is needed to calculate aggregation over contigs below.
        # This is the step where the genome is split into its contigs for parallel execution.
        # The get_fai() function uses a snakemake checkpoint to make sure that the fai is
        # produced before we use it here to get its content.
        ref=get_fai,
        contig_groups=contigs_groups_input,
        gvcfs=combined_contig_gvcfs,
        dones=combined_contig_done,
    output:
        # vcf="calling/genotyped-all.vcf.gz",
        # tbi="calling/genotyped/all.vcf.gz.tbi",
        # lst="calling/genotyped/all.txt",
        # done="calling/genotyped-all.done"
        vcf=(
            temp("calling/merged/merged-all.vcf.gz")
            if (config["settings"].get("contig-group-size"))
            else "calling/genotyped-all.vcf.gz"
        ),
        tbi=(
            temp("calling/merged/merged-all.vcf.gz.tbi")
            if (config["settings"].get("contig-group-size"))
            else "calling/genotyped-all.vcf.gz.tbi"
        ),
        lst=(
            temp("calling/merged/merged-all.txt")
            if (config["settings"].get("contig-group-size"))
            else "calling/genotyped-all.txt"
        ),
        done=(
            touch("calling/merged/merged-all.vcf.gz.done")
            if (config["settings"].get("contig-group-size"))
            else touch("calling/genotyped-all.vcf.gz.done")
        ),
    params:
        # Use a list of the chromosomes in the same order as the fai, for bcftools to sort the output.
        regionorder=combined_contig_order,
    log:
        "logs/calling/bcftools-combine-all.log",
    benchmark:
        "benchmarks/calling/bcftools-combine-all.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        # Create a list of the input files. We use the order of gvfcs as provided by the fai file
        # by default, which works if we call per chromosome, so that bcftool has an easier job
        # with the data, but gives them out of order if we use  contig groups insteads.
        # Hence, we further enforce the order that we want by specifing the chromosomes as regions,
        # see https://github.com/samtools/bcftools/issues/268#issuecomment-105772010,
        # by providing a list of regions in the order that we want.
        "echo {input.gvcfs} | sed 's/ /\\n/g' > {output.lst} ; "
        "("
        "bcftools concat "
        '   --file-list "{output.lst}" --regions "{params.regionorder}" --allow-overlaps '
        "   --output-type z -o {output.vcf} ; "
        "bcftools index --tbi {output.vcf}"
        ") &> {log}"
