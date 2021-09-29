# =================================================================================================
#     Common Helper Functions
# =================================================================================================

def get_fai(wildcards):
    # Stop at the snakemake checkpoint first to ensure that the fai file is available.
    return checkpoints.samtools_faidx.get().output[0]
    # return config["data"]["reference"]["genome"] + ".fai"

# Contigs in reference genome.
def get_contigs( fai ):
    small_contig_thresh = config["settings"].get("small-contigs-threshold", 0)
    if small_contig_thresh > 0:
        contigs = []
        length_sum = 0
        contig0 = ""
        with open(fai, "r") as f:
            for line in f:
                contig, length_str = line.split("\t")[:2]
                contig = contig.strip()
                length = int(length_str.strip())

                if not contig0:
                    contig0 = contig

                length_sum += length

                if length_sum >= small_contig_thresh:
                    contigs.append(f"{contig0}---{contig}")
                    contig0 = ""
                    length_sum = 0

        if contig0:
            contigs.append(f"{contig0}---{contig}")

        return contigs

    return pd.read_csv( fai, sep='\t', header=None, usecols=[0], squeeze=True, dtype=str )

# =================================================================================================
#     Restrict Regions & Small Contigs
# =================================================================================================

# Interset the restict regions file with a given contig (chromosome), so that we can use the
# resulting bed file for parallelization over contigs.
if "restrict-regions" in config["settings"]:
    rule compose_regions:
        input:
            config["settings"].get("restrict-regions")
        output:
            "called/{contig}.regions.bed"
        log:
            "logs/bedextract/{contig}.regions.log"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"

    # Rule is not submitted as a job to the cluster.
    localrules: compose_regions

if config["settings"].get("small-contigs-threshold", 0) > 0:
    rule group_contigs:
        input:
            ref=get_fai
        output:
            "contig-groups/{contig}.list"
        log:
            "logs/contig-groups/{contig}.log"
        run:
            contig0, contig1 = wildcards.contig.split("---")
            contigs = []
            with open(input.ref, "r") as f:
                for line in f:
                    contig = line.split("\t")[0].strip()
                    if contigs or contig == contig0:
                        contigs.append(contig)
                    if contig == contig1:
                        break

            with open(output[0], "w") as f:
                f.writelines(f"{c}\n" for c in contigs)

    # Rule is not submitted as a job to the cluster.
    localrules: group_contigs

    # Conflicts of interest:
    if config["settings"].get("restrict-regions"):
        raise Exception(
            "Cannot combine settings small-contigs-threshold > 0 with restrict-regions "
            "at the moment, as we have not implemented this yet. "
            "If you need this combination of settings, please submit an issue to "
            "https://github.com/lczech/grenepipe/issues and we will see what we can do."
        )

    if config["settings"]["calling-tool"] != "haplotypecaller":
        raise Exception(
            "Can only use setting small-contigs-threshold with calling-tool haplotypecaller "
            "at the moment, as we have not implemented this for other calling tools yet. "
            "If you need this combination of settings, please submit an issue to "
            "https://github.com/lczech/grenepipe/issues and we will see what we can do."
        )

# =================================================================================================
#     Variant Calling
# =================================================================================================

# Switch to the chosen caller
if config["settings"]["calling-tool"] == "haplotypecaller":

    # Use `GATK HaplotypeCaller`
    include: "calling-haplotypecaller.smk"

elif config["settings"]["calling-tool"] == "bcftools":

    # Use `bcftools call`
    include: "calling-bcftools.smk"

elif config["settings"]["calling-tool"] == "freebayes":

    # Use `freebayes`
    include: "calling-freebayes.smk"

else:
    raise Exception("Unknown calling-tool: " + config["settings"]["calling-tool"])
