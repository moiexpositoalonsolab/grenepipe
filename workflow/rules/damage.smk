import platform

# =================================================================================================
#     mapDamage
# =================================================================================================


rule mapdamage:
    input:
        # Get either the normal mapping output, or, if additional filtering via `samtools view`
        # is set in the config settings: filter-mapped-reads, use the filtered output instead.
        get_mapped_reads,
        # "mapping/sorted/{sample}-{unit}.bam"
        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
    output:
        "damage/mapdamage/{sample}/Runtime_log.txt",
    params:
        index=config["data"]["reference-genome"],
        extra=config["params"]["mapdamage"]["extra"],
        outdir="damage/mapdamage/{sample}",
    log:
        "logs/damage/mapdamage/{sample}.log",
    conda:
        # We have two different env yaml files, depending on the platform.
        # This is because on Linux, a particular lib might be missing that is needed for mapdamage,
        # but that lib does not exist on conda for MacOS, so we want to skip it there.
        (
            "../envs/mapdamage-macos.yaml"
            if platform.system() == "Darwin"
            else "../envs/mapdamage-linux.yaml"
        )
    shell:
        "mapDamage -i {input[0]} -r {params.index} -d {params.outdir} {params.extra} > {log} 2>&1"


rule mapdamage_collect:
    input:
        expand("damage/mapdamage/{sample}/Runtime_log.txt", sample=config["global"]["sample-names"]),
    output:
        touch("damage/mapdamage/mapdamage.done"),


localrules:
    mapdamage_collect,


# =================================================================================================
#     DamageProfiler
# =================================================================================================


rule damageprofiler:
    input:
        # Get either the normal mapping output, or, if additional filtering via `samtools view`
        # is set in the config settings: filter-mapped-reads, use the filtered output instead.
        get_mapped_reads,
        # "mapping/sorted/{sample}-{unit}.bam"
        # Get the reference genome and its indices. Not sure if the indices are needed
        # for this particular rule, but doesn't hurt to include them as an input anyway.
        ref=config["data"]["reference-genome"],
        refidcs=expand(
            config["data"]["reference-genome"] + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa", "fai"],
        ),
    output:
        "damage/damageprofiler/{sample}/DamageProfiler.log",
    params:
        index=config["data"]["reference-genome"],
        extra=config["params"]["damageprofiler"]["extra"],
        outdir="damage/damageprofiler/{sample}",
    log:
        "logs/damage/damageprofiler/{sample}.log",
    conda:
        "../envs/damageprofiler.yaml"
    shell:
        "damageprofiler -i {input[0]} -r {params.index} -o {params.outdir} {params.extra} > {log} 2>&1"
        # "java -jar DamageProfiler-0.5.0.jar "


rule damageprofiler_collect:
    input:
        expand(
            "damage/damageprofiler/{sample}/DamageProfiler.log",
            sample=config["global"]["sample-names"],
        ),
    output:
        touch("damage/damageprofiler/damageprofiler.done"),


localrules:
    damageprofiler_collect,
