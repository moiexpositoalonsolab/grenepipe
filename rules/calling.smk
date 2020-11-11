# =================================================================================================
#     Restrict Regions
# =================================================================================================

# Interset the restict regions file with a given contig (chromosome), so that we can sue the
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

# =================================================================================================
#     Common Helper Functions
# =================================================================================================

def get_fai():
    return config["data"]["reference"]["genome"] + ".fai"

# contigs in reference genome
def get_contigs():
    return pd.read_csv(get_fai(), sep='\t', header=None, usecols=[0], squeeze=True, dtype=str)

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
