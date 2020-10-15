# =================================================================================================
#     mapDamage
# =================================================================================================

rule mapdamage:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        "mapdamage/{sample}-{unit}/Runtime_log.txt"
    params:
        index=config["data"]["reference"]["genome"],
        extra=config["params"]["mapdamage"]["extra"],
        outdir="mapdamage/{sample}-{unit}"
    log:
        "logs/mapdamage/{sample}-{unit}.log"
    conda:
        "../envs/mapdamage.yaml"
    shell:
        "mapDamage -i {input} -r {params.index} -d {params.outdir} {params.extra} > {log} 2>&1"

# =================================================================================================
#     DamageProfiler
# =================================================================================================

rule damageprofiler:
    input:
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        "damageprofiler/{sample}-{unit}/DamageProfiler.log"
    params:
        index=config["data"]["reference"]["genome"],
        extra=config["params"]["damageprofiler"]["extra"],
        outdir="damageprofiler/{sample}-{unit}",

        # Unfortunately, the tool creates another output folder below the specified one,
        # contraticting its documentation. Grrrrr. So, we have to move the files ourselves...
        bogus_outdir="damageprofiler/{sample}-{unit}/{sample}-{unit}.sorted"
    log:
        "logs/damageprofiler/{sample}-{unit}.log"
    conda:
        "../envs/damageprofiler.yaml"
    shell:
        "damageprofiler -i {input} -r {params.index} -o {params.outdir} "
        "{params.extra} > {log} 2>&1 ;"
        "mv {params.bogus_outdir}/* {params.outdir}; "
        "rm -r {params.bogus_outdir}"
        # "java -jar DamageProfiler-0.5.0.jar "
