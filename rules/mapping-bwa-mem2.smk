# =================================================================================================
#     Read Mapping
# =================================================================================================

rule map_reads:
    input:
        reads=get_trimmed_reads,
        ref=config["data"]["reference"]["genome"],
        refidcs=expand(
            config["data"]["reference"]["genome"] + ".{ext}",
            ext=[ "0123", "amb", "ann", "bwt.2bit.64", "pac" ]
        )
    output:
        "mapped/{sample}-{unit}.sorted.bam"
    params:
        index=config["data"]["reference"]["genome"],

        # We need the read group tags, including `ID` and `SM`, as downstream tools use these.
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}' " + config["params"]["bwamem2"]["extra"],

        # Sort as we need it.
        sort="samtools",
        sort_order="coordinate",
        sort_extra=config["params"]["bwamem2"]["extra-sort"]
    group:
        "mapping"
    log:
        "logs/bwa-mem2/{sample}-{unit}.log"
    benchmark:
        "benchmarks/bwa-mem2/{sample}-{unit}.bench.log"
    threads:
        config["params"]["bwamem2"]["threads"]
    conda:
        # As always, we need our own env here that overwrites the python/pandas/numpy stack
        # to make sure that we do not run into a version conflict.
        "../envs/bwa-mem2.yaml"
    wrapper:
        "0.78.0/bio/bwa-mem2/mem"
