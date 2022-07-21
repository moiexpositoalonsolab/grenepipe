# =================================================================================================
#     Mark Duplicates
# =================================================================================================

# Used in `mapping.smk`

# The dedup documentation says that -o is an output file, but the program then complains
# that it needs a directory instead... and then litters the output directory with all kinds
# of files that we do not want. So, let's use a shadow directory, and let's pray that the file naming
# of the `*rmdup.bam` file is consistent (hard to tell, as this is not documented in dedup)...
# --> Update: Not using a shadow diretory any more, still seems to work fine.
# Furthermore, dedup does not sort, so we have to do this ourselves (Again! We already sorted
# after mapping... How does that program mess up the existing order?!).
# Luckily, this combination of unfortunate usages (weird output files, and additional need to sort)
# can be leveraged here by declaring the unsorted, weirdly named dedup output file as a temp file,
# and then feed that into samtools for sorting, after which the dedup file will be deleted again.
# At least a little bit of comfort!
# Lastly, we also keep the json file for reporting with multiqc.
rule mark_duplicates:
    input:
        # Get either the normal mapping output, or, if additional filtering via `samtools view`
        # is set in the config settings: filter-mapped-reads, use the filtered output instead.
        get_mapped_reads
        # "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=temp("dedup/{sample}-{unit}_rmdup.bam"),
        metrics="dedup/{sample}-{unit}.dedup.json"
    log:
        "logs/dedup/{sample}-{unit}.log"
    benchmark:
        "benchmarks/dedup/{sample}-{unit}.bench.log"
    params:
        extra=config["params"]["dedup"]["extra"],
        out_dir="dedup"
    conda:
        "../envs/dedup.yaml"
    group:
        "mapping_extra"
    shell:
        # Dedup bases its output names on this, so we need some more trickery to make it
        # output the file names that we want.
        "dedup -i {input} -o {params.out_dir} {params.extra} > {log} 2>&1 ; "
        "mv"
        "    dedup/{wildcards.sample}-{wildcards.unit}." + get_mapped_read_infix() + "_rmdup.bam"
        "    dedup/{wildcards.sample}-{wildcards.unit}_rmdup.bam ; "
        "mv"
        "    dedup/{wildcards.sample}-{wildcards.unit}." + get_mapped_read_infix() + ".dedup.json"
        "    dedup/{wildcards.sample}-{wildcards.unit}.dedup.json"

rule sort_reads_dedup:
    input:
        "dedup/{sample}-{unit}_rmdup.bam"
    output:
        (
            "dedup/{sample}-{unit}.bam"
            if config["settings"]["keep-intermediate"]["mapping"]
            else temp("dedup/{sample}-{unit}.bam")
        )
    params:
        extra=config["params"]["samtools"]["sort"],
        tmp_dir=config["params"]["samtools"]["temp-dir"]
    threads:  # Samtools takes additional threads through its option -@
        1     # This value - 1 will be sent to -@. Weird flex, but okay.
    log:
        "logs/samtools/sort/{sample}-{unit}-dedup.log"
    group:
        "mapping_extra"
    wrapper:
        "0.80.0/bio/samtools/sort"
