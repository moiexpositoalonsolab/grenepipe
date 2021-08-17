# =================================================================================================
#     Variant Calling
# =================================================================================================

if config["data"]["reference"].get("known-variants"):
    raise Exception("Calling tool 'bcftools' cannot be used with the option 'known-variants'")

def get_mpileup_params(wildcards, input):
    # Start with the user-specified params from the config file
    params = config["params"]["bcftools"]["mpileup"]

    # input.regions is either the restrict-regions file, or empty.
    # If given, we need to set a param of mpileup to use the file.
    # If not given, we set a param that uses the contig string instead.
    if input.regions:
        params += " --regions-file {}".format(input.regions)
    else:
        params += " --regions {}".format(wildcards.contig)

    # Lastly, we need some extra annotations, so that our downstream stats rules can do their job.
    # That was a hard debugging session to figure this one out... If only bioinformatics tools
    # had better error messages!
    params += " --annotate FORMAT/AD,FORMAT/DP"
    return params

rule call_variants:
    input:
        # Need the ref genome, as well as its indices.
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),

        # Get the bam and bai files for all files. If this is too slow, we need to split,
        # similar to what our GATK HaplotypeCaller rule does.
        samples=get_all_bams(),
        indices=get_all_bais(),

        # If we use restricted regions, set them here. If not, empty, which will propagate to the
        # get_mpileup_params function as well
        regions="called/{contig}.regions.bed" if config["settings"].get("restrict-regions") else []
    output:
        vcf=protected("called/{contig}.vcf.gz")
    params:
        # Optional parameters for bcftools mpileup (except -g, -f).
        mpileup=get_mpileup_params,

        # Optional parameters for bcftools call (except -v, -o, -m).
        call=config["params"]["bcftools"]["call"]
    log:
        "logs/bcftools/call-{contig}.log"
    benchmark:
        "benchmarks/bcftools/call-{contig}.bench.log"
    conda:
        "../envs/bcftools.yaml"
    threads:
        config["params"]["bcftools"]["threads"]
    shell:
        "(bcftools mpileup {params.mpileup} --fasta-ref {input.ref} --output-type u {input.samples} | "
        "bcftools call --threads {threads} {params.call} --output-type z -o {output[0]} -v -) 2> {log}"

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
# call it mergeing here...

# Need an input function to work with the fai checkpoint
def merge_variants_vcfs_input(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("called/{contig}.vcf.gz", contig=get_contigs( fai ))

# Need index files for some of the downstream tools.
# Rational for the fai: see merge_variants_vcfs_input()
def merge_variants_tbis_input(wildcards):
    fai = checkpoints.samtools_faidx.get().output[0]
    return expand("called/{contig}.vcf.gz.tbi", contig=get_contigs( fai ))

# bcftools does not automatically create vcf index files, so we need a rule for that...
# ... but picard does! So, if we continue using the above picard mergevcfs tool, we don't need tabix
# ...... buuuuut vcflib does not! So, we reactivate it again!
rule vcf_index:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf"
    log:
        "logs/tabix/{prefix}.log"
    wrapper:
        "0.55.1/bio/tabix"

# Cannot use bcftools concat, as it does not except vcf/bcf files without any calls in them.
# rule merge_variants:
#     input:
#         calls=lambda w: expand("called/{contig}.bcf", contig=get_contigs())
#     output:
#         "genotyped/all.vcf.gz"
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
#         # vcfs=lambda w: expand("called/{contig}.vcf.gz", contig=get_contigs())
#         vcfs=merge_variants_vcfs_input
#     output:
#         vcf="genotyped/all.vcf.gz"
#     log:
#         "logs/picard/merge-genotyped.log"
#     benchmark:
#         "benchmarks/picard/merge-genotyped.bench.log"
#     conda:
#         "../envs/picard.yaml"
#     wrapper:
#         "0.51.3/bio/picard/mergevcfs"

# Let's try vcflib instead, hoping that it can handle both of the above use cases.
# See above for comments on the input files.
# This is not the safest option, as vcflib does little to no checks for file integrity,
# for example that the headers and sample columsn actually are good to be merged.
# But as we are using this in our pipeline, where we know that the data is good, this should be fine.
rule merge_variants:
    input:
        ref=get_fai,
        vcfs=merge_variants_vcfs_input,
        tbis=merge_variants_tbis_input
    output:
        vcf="genotyped/all.vcf.gz"
    log:
        "logs/vcflib/merge-genotyped.log"
    benchmark:
        "benchmarks/vcflib/merge-genotyped.bench.log"
    conda:
        "../envs/vcflib.yaml"
    shell:
        "vcfcat {input.vcfs} | bgzip --stdout > {output.vcf} && tabix -p vcf {output.vcf}"
