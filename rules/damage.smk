# =================================================================================================
#     mapDamage
# =================================================================================================

rule mapdamage:
    input:
        # Get either the normal mapping output, or, if additional filtering via `samtools view`
        # is set in the config settings: filter-mapped-reads, use the filtered output instead.
        get_mapped_reads,
        # "mapped/{sample}-{unit}.sorted.bam"

        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),
    output:
        "mapdamage/{sample}-{unit}/Runtime_log.txt"
    params:
        index=config["data"]["reference-genome"],
        extra=config["params"]["mapdamage"]["extra"],
        outdir="mapdamage/{sample}-{unit}"
    log:
        "logs/mapdamage/{sample}-{unit}.log"
    conda:
        "../envs/mapdamage.yaml"
    shell:
        "mapDamage -i {input[0]} -r {params.index} -d {params.outdir} {params.extra} > {log} 2>&1"

rule mapdamage_collect:
    input:
        expand(
            "mapdamage/{u.sample}-{u.unit}/Runtime_log.txt",
            u=config["global"]["samples"].itertuples()
        )
    output:
        touch("mapdamage/mapdamage.done")

localrules: mapdamage_collect

# =================================================================================================
#     DamageProfiler
# =================================================================================================

rule damageprofiler:
    input:
        # Get either the normal mapping output, or, if additional filtering via `samtools view`
        # is set in the config settings: filter-mapped-reads, use the filtered output instead.
        get_mapped_reads,
        # "mapped/{sample}-{unit}.sorted.bam"

        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),
    output:
        "damageprofiler/{sample}-{unit}/DamageProfiler.log"
    params:
        index=config["data"]["reference-genome"],
        extra=config["params"]["damageprofiler"]["extra"],
        outdir="damageprofiler/{sample}-{unit}"
    log:
        "logs/damageprofiler/{sample}-{unit}.log"
    conda:
        "../envs/damageprofiler.yaml"
    shell:
        "damageprofiler -i {input[0]} -r {params.index} -o {params.outdir} {params.extra} > {log} 2>&1"
        # "java -jar DamageProfiler-0.5.0.jar "

rule damageprofiler_collect:
    input:
        expand(
            "damageprofiler/{u.sample}-{u.unit}/DamageProfiler.log",
            u=config["global"]["samples"].itertuples()
        )
    output:
        touch("damageprofiler/damageprofiler.done")

localrules: damageprofiler_collect
