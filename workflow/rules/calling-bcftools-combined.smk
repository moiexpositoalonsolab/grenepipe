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
        # Get the bam and bai files for all files. If this is too slow, we need to split,
        # similar to what our GATK HaplotypeCaller rule does.
        samples=get_all_bams(),
        indices=get_all_bais(),
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
        vcf=(
            "calling/called/{contig}.vcf.gz"
            if config["settings"]["keep-intermediate"]["calling"]
            else temp("calling/called/{contig}.vcf.gz")
        ),
        # vcf=protected("calling/called/{contig}.vcf.gz")
        done=touch("calling/called/{contig}.done"),
    params:
        # Optional parameters for bcftools mpileup (except -g, -f).
        mpileup=get_mpileup_params,
        # Optional parameters for bcftools call (except -v, -o, -m).
        call=config["params"]["bcftools"]["call"],
    log:
        "logs/calling/bcftools/call-{contig}.log",
    benchmark:
        "benchmarks/calling/called/bcftools/call-{contig}.log"
    conda:
        "../envs/bcftools.yaml"
    threads: config["params"]["bcftools"]["threads"]
    shell:
        # We run an additional norm step after the calling, in order to clear any calls that
        # contain ambiguity chars. As always, no idea how they end up in the vcf in the first place.
        # See https://github.com/samtools/bcftools/issues/473 for details.
        "("
        "bcftools mpileup {params.mpileup} --fasta-ref {input.ref} --output-type u {input.samples} | "
        "bcftools call --threads {threads} {params.call} --variants-only --output-type u - | "
        "bcftools norm --check-ref ws --fasta-ref {input.ref} --output-type z --output {output.vcf} - "
        ") 2> {log}"


# Cannot use the wrapper, as we need `bcftools mpileup` instead of `samtools mpileup` (as used
# by the wrapper), in order to have the `--annotate` option, which we need so that our
# stats rules can work. Man... what a mess! Again!
# wrapper:
#     "0.55.1/bio/bcftools/call"

# Alternative to parallelization over contigs would be to parallelize over equally sized chunks of
# the reference genome, see e.g., http://dmnfarrell.github.io/bioinformatics/mpileup-parallel

# =================================================================================================
#     Combining Calls
# =================================================================================================

# Technically, this step should be a "concatenation" (combine multiple vcfs of the same individuals,
# but different chromosomes/contigs), not a "merging" (combine multiple vcfs of the same contig,
# but with different individuals)... But tools tend to name these things differently, so let's
# call it merging here...


# Need an input function to work with the fai checkpoint
def merge_variants_vcfs_input(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("calling/called/{contig}.vcf.gz", contig=get_contigs(fai))


# Need index files for some of the downstream tools.
# Rational for the fai: see merge_variants_vcfs_input()
def merge_variants_tbis_input(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("calling/called/{contig}.vcf.gz.tbi", contig=get_contigs(fai))


# bcftools does not automatically create vcf index files, so we need a rule for that...
# ... but picard does! So, if we continue using the below picard mergevcfs tool, we don't need tabix
# ...... buuuuut vcflib does not! So, we reactivate it again!
rule called_vcf_index:
    input:
        "calling/called/{contig}.vcf.gz",
    output:
        "calling/called/{contig}.vcf.gz.tbi",
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf",
    log:
        "logs/calling/tabix/{contig}.log",
    conda:
        "../envs/tabix.yaml"
    wrapper:
        "0.55.1/bio/tabix"


# Cannot use bcftools concat, as it does not except vcf/bcf files without any calls in them.
# rule merge_variants:
#     input:
#         calls=lambda w: expand("calling/called/{contig}.bcf", contig=get_contigs())
#     output:
#         "calling/genotyped-all.vcf.gz"
#     log:
#         "logs/bcftools/concat.log"
#     params:
#         # optional parameters for bcftools concat (except -o)
#         "--output-type z " + config["params"]["bcftools"]["call"]
#     wrapper:
#         "0.55.1/bio/bcftools/concat"

# Can also not use picard, as it does not handle the ambiguity codes in variant calls
# that `bcftools call` can produce...
# rule merge_variants:
#     input:
#         # fai is needed to calculate aggregation over contigs below.
#         # This is the step where the genome is split into its contigs for parallel execution.
#         # The get_fai() function uses a snakemake checkpoint to make sure that the fai is
#         # produced before we use it here to get its content.
#         ref=get_fai,
#
#         # The wrapper expects input to be called `vcfs`, but we can use `vcf.gz` as well.
#         # vcfs=lambda w: expand("calling/called/{contig}.vcf.gz", contig=get_contigs())
#         vcfs=merge_variants_vcfs_input
#     output:
#         vcf="calling/genotyped-all.vcf.gz"
#     log:
#         "logs/picard/merge-genotyped.log"
#     benchmark:
#         "benchmarks/picard/merge-genotyped.log"
#     conda:
#         "../envs/picard.yaml"
#     wrapper:
#         "0.51.3/bio/picard/mergevcfs"


# Let's try vcflib instead, hoping that it can handle both of the above use cases.
# See above for comments on the input files.
# This is not the safest option, as vcflib does little to no checks for file integrity,
# for example that the headers and sample columns actually are good to be merged.
# But as we are using this in our pipeline, where we know that the data is good, this should be fine.
rule merge_variants:
    input:
        ref=get_fai,
        contig_groups=contigs_groups_input,
        vcfs=merge_variants_vcfs_input,
        tbis=merge_variants_tbis_input,
    output:
        # With small contigs, we also need sorting, see below.
        # Unfortunately, we cannot pipe here, as Picard fails with that, so temp file it is...
        # If we do not use small contigs, we directly output the final file.
        vcf=(
            temp("calling/merged/merged-all.vcf.gz")
            if (config["settings"].get("contig-group-size"))
            else "calling/genotyped-all.vcf.gz"
        ),
        done=(
            touch("calling/merged/merged-all.done")
            if (config["settings"].get("contig-group-size"))
            else touch("calling/genotyped-all.done")
        ),
    log:
        "logs/calling/vcflib/merge-genotyped.log",
    benchmark:
        "benchmarks/calling/genotyped/vcflib/merge-genotyped.log"
    conda:
        "../envs/vcflib.yaml"
    shell:
        "vcfcat {input.vcfs} | bgzip --stdout > {output.vcf} && tabix -p vcf {output.vcf}"
