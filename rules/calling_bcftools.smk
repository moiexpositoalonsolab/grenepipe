# =================================================================================================
#     Helper Functions
# =================================================================================================

from itertools import chain

# Return the bam file(s) for a given sample
def get_sample_bams(sample):
    return expand(get_mapping_result(), sample=sample, unit=samples.loc[sample].unit)

# Return the bai file(s) for a given sample
def get_sample_bais(sample):
    return expand(get_mapping_result(True), sample=sample, unit=samples.loc[sample].unit)

# Return the bam file(s) for all samples
def get_all_bams():
    return list(chain.from_iterable( [ get_sample_bams(sample) for sample in sample_names ] ))

# Return the bai file(s) for all samples
def get_all_bais():
    return list(chain.from_iterable( [ get_sample_bais(sample) for sample in sample_names ] ))

# =================================================================================================
#     Variant Calling
# =================================================================================================

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
        ref=config["data"]["reference"]["genome"],

        # Get the bam and bai files for the given sample only.
        samples=get_all_bams(),
        indices=get_all_bais(),

        # If we use restricted regions, set them here. If not, empty, which will propagate to the
        # get_mpileup_params function as well
        regions=config["rundir"] + "called/{contig}.regions.bed" if config["settings"].get("restrict-regions") else []
    output:
        vcf=protected(config["rundir"] + "called/{contig}.vcf.gz")
    params:
        # Optional parameters for bcftools mpileup (except -g, -f).
        mpileup=get_mpileup_params,

        # Optional parameters for bcftools call (except -v, -o, -m).
        call=config["params"]["bcftools"]["call"]
    log:
        config["rundir"] + "logs/bcftools/call-{contig}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools mpileup {params.mpileup} --fasta-ref {input.ref} --output-type u {input.samples} | "
        "bcftools call -m {params.call} --output-type z -o {output[0]} -v -) 2> {log}"

    # Cannot use the wrapper, as we need `bcftools mpileup` instead of `samtools mpileup` (as used
    # by the wrapper), in order to have the `--annotate` option, which we need so that our
    # stats rules can work.
    # wrapper:
    #     "0.55.1/bio/bcftools/call"

# Alternative to parallelization over contigs would be to parallelize over equally sized chunks of
# the reference genome, see e.g., http://dmnfarrell.github.io/bioinformatics/mpileup-parallel

# =================================================================================================
#     Combining Calls
# =================================================================================================

# Cannot use bcftools concat, as it does not except vcf/bcf files without any calls in them.
# rule merge_variants:
#     input:
#         calls=lambda w: expand(config["rundir"] + "called/{contig}.bcf", contig=get_contigs())
#     output:
#         config["rundir"] + "genotyped/all.vcf.gz"
#     log:
#         config["rundir"] + "logs/bcftools/concat.log"
#     params:
#         # optional parameters for bcftools concat (except -o)
#         "--output-type z " + config["params"]["bcftools"]["call"]
#     wrapper:
#         "0.55.1/bio/bcftools/concat"

rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below

        # The wrapper expects input to be called `vcfs`, but we can use `vcf.gz` as well.
        vcfs=lambda w: expand(config["rundir"] + "called/{contig}.vcf.gz", contig=get_contigs())
    output:
        vcf=config["rundir"] + "genotyped/all.vcf.gz"
    log:
        config["rundir"] + "logs/picard/merge-genotyped.log"
    wrapper:
        "0.51.3/bio/picard/mergevcfs"

# bcftools does not automatically create vcf index files, so we need a rule for that...
# ... but picard does! So, if we continue using the above picard mergevcfs tool, we don't need tabix
# rule vcf_index:
#     input:
#         "{prefix}.vcf.gz"
#     output:
#         "{prefix}.vcf.gz.tbi"
#     params:
#         # pass arguments to tabix (e.g. index a vcf)
#         "-p vcf"
#     log:
#             config["rundir"] + "logs/tabix/{prefix}.log"
#     wrapper:
#         "0.55.1/bio/tabix"
