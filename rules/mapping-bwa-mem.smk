# =================================================================================================
#     Read Mapping
# =================================================================================================

# `samtools sort` creates tmp files that are not cleaned up when the cluster job runs out
# of time, but which cause samtools to immediately terminate if called again, meaning
# that we cannot run it again with more wall time. Hence, we need to set a specific tmp dir,
# which we then let snakemake clean up for us by marking it as temp. The `samtools sort` manual
# states that when given an existing directory to the `-T` option to specify the temp output,
# it uses tmp file names that are unique per run, so we can re-use the same dir for each attempt
# at least, without collision. But we need to make individual dirs per sample/unit, so that the
# automatic deletion can do its job on a per-sample/unit basis.
rule mk_samtools_sort_tmp_dir:
    output:
        tmpdir=temp(directory( config["rundir"] + "mapped/samtools-sort-{sample}-{unit}/" ))
    shell:
        "mkdir -p {output.tmpdir}"

localrules: mk_samtools_sort_tmp_dir

rule map_reads:
    input:
        reads=get_trimmed_reads,
        tmpdir=config["rundir"] + "mapped/samtools-sort-{sample}-{unit}/"
    output:
        config["rundir"] + "mapped/{sample}-{unit}.sorted.bam"
    log:
        config["rundir"] + "logs/bwa-mem/{sample}-{unit}.log"
    benchmark:
        config["rundir"] + "benchmarks/bwa-mem/{sample}-{unit}.bench.log"
    params:
        index=config["data"]["reference"]["genome"],

        # We need the read group tags, including `ID` and `SM`, as downstream tools use these.
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        # TODO add LD field as well for the unit?! http://www.htslib.org/workflow/

        sort="samtools",         # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate", # Can be 'queryname' or 'coordinate'.
        # Extra args for samtools/picard.
        sort_extra=" -T {input.tmpdir} "
    threads:
        config["params"]["bwamem"]["threads"]
    resources:
        # Increase time limit in factors of 2h, if the job fails due to time limit.
        time = lambda wildcards, attempt: int(120 * int(attempt))
    shadow: "full"
    wrapper:
        "0.51.3/bio/bwa/mem"
