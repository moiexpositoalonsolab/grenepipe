# =================================================================================================
#     Variant Calling
# =================================================================================================

rule call_variants:
    input:
        # Need the ref genome, as well as its indices.
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),

        # Get the sample data.
        samples=get_sample_bams_wildcard,
        indices=get_sample_bais_wildcard,

        # If we use restricted regions, set them here. If not, empty, which will propagate to the
        # get_mpileup_params function as well. Same for small contig groups.
        # regions="called/{contig}.regions.bed" if config["settings"].get("restrict-regions") else []
        regions="called/{contig}.regions.bed" if (
            config["settings"].get("restrict-regions")
        ) else (
            "contig-groups/{contig}.bed" if (
                config["settings"].get("contig-group-size")
            ) else []
        )
    output:
        gvcf="called/{sample}.{contig}.g.vcf.gz",
        gtbi="called/{sample}.{contig}.g.vcf.gz.tbi"
    params:
        # Optional parameters for bcftools mpileup (except -g, -f).
        mpileup=get_mpileup_params,

        # Optional parameters for bcftools call (except -v, -o, -m).
        call=config["params"]["bcftools"]["call"]
    log:
        "logs/bcftools/call-{sample}.{contig}.log"
    benchmark:
        "benchmarks/bcftools/call-{sample}.{contig}.bench.log"
    conda:
        "../envs/bcftools.yaml"
    threads:
        config["params"]["bcftools"]["threads"]
    shell:
        "("
        "bcftools mpileup {params.mpileup} --fasta-ref {input.ref} --output-type u {input.samples} "
        "-a 'INFO/AD' -a 'FORMAT/AD' -a 'FORMAT/DP' --gvcf 5,10,25,50 | "
        "bcftools call --threads {threads} {params.call} --gvcf 5,10,25,50 "
        "--output-type z -o {output.gvcf} ; "
        "bcftools index --tbi {output.gvcf}"
        ") &> {log}"

# =================================================================================================
#     Combining Contigs
# =================================================================================================

rule combine_contig:
    input:
        # Need the ref genome, as well as its indices.
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),

        # Get the sample data, including indices, which are produced above.
        gvcfs=expand(
            "called/{sample}.{{contig}}.g.vcf.gz",
            sample=config["global"]["sample-names"]
        ),
        indices=expand(
            "called/{sample}.{{contig}}.g.vcf.gz.tbi",
            sample=config["global"]["sample-names"]
        )
    output:
        gvcf="called/all.{contig}.g.vcf.gz",
        gtbi="called/all.{contig}.g.vcf.gz.tbi",
        gvcflist="called/all.{contig}.g.txt"
    log:
        "logs/bcftools/combine-contig-{contig}.log"
    benchmark:
        "benchmarks/bcftools/combine-contig-{contig}.bench.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        # Sort the input gvcfs into a list file, so that we always have the same order,
        # in the hope that this works with large numbers of files without exceeding the ulimit
        # for number of files that can be kept open. Needs to be tested on a large dataset though...
        "echo {input.gvcfs} | sed 's/ /\\n/g' | sort > {output.gvcflist} ;"
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

# Need an input function to work with the fai checkpoint
def combined_contig_gvcfs(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("called/all.{contig}.g.vcf.gz", contig=get_contigs( fai ))

# We also need a comma-separated list of the contigs, so that bcftools can output
# the concatenated entries in the correct order as given in the fai.
# For this, we use the same technique of using the fai checkpoint as before.
def combined_contig_order(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    contig_list = []
    with open(fai, "r") as f:
        for line in f:
            contig = line.split("\t")[0].strip()
            contig_list.append( contig )
    return ",".join(contig_list)

rule combine_all:
    input:
        # fai is needed to calculate aggregation over contigs below.
        # This is the step where the genome is split into its contigs for parallel execution.
        # The get_fai() function uses a snakemake checkpoint to make sure that the fai is
        # produced before we use it here to get its content.
        ref=get_fai,
        gvcfs=combined_contig_gvcfs
    output:
        vcf="genotyped/all.vcf.gz",
        tbi="genotyped/all.vcf.gz.tbi",
        lst="genotyped/all.txt"
    params:
        regionorder=combined_contig_order
    log:
        "logs/bcftools/combine-all.log"
    benchmark:
        "benchmarks/bcftools/combine-all.bench.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        # Use the order of gvfcs as provided by the fai file by default, which works if we call
        # per chromosome, but gives them out of order if we use contig groups insteads.
        "echo {input.gvcfs} | sed 's/ /\\n/g' > {output.lst} ; "
        # Hence enforce the order by specifing them as regions, see
        # https://github.com/samtools/bcftools/issues/268#issuecomment-105772010
        # by providing a list of regions in the order that we want.
        "("
        "bcftools concat --file-list \"{output.lst}\" --regions \"{params.regionorder}\" --allow-overlaps "
        "--output-type z -o {output.vcf} ;"
        "bcftools index --tbi {output.vcf}"
        ") &> {log}"
