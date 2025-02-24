# =================================================================================================
#     Grouping of (Small) Contigs
# =================================================================================================


# Get the list of chromosome names that are present in the fai file,
# and their length, with a length filter if needed.
def read_contigs_from_fai(fai, min_contig_size=0):
    # return list(pd.read_csv(fai, sep="\t", header=None, usecols=[0], dtype=str).squeeze())

    # Read fai to get all contigs and their sizes.
    contig_list = []
    with open(fai, "r") as f:
        for line in f:
            contig, length_str = line.split("\t")[:2]
            contig = contig.strip()
            length = int(length_str.strip())
            if length < min_contig_size:
                continue
            contig_list.append((contig, length))
    return contig_list


# If we want to combine contigs into groups, use the rules and functions for this.
if config["settings"].get("contig-group-size", 0) > 0:

    include: "calling-contig-groups.smk"


# Dummy definition of the above rule for when we are not using contig groups.
# We cannot do this as an else branch of the above, due to an issue in the
# formatting of snakemfmt: https://github.com/snakemake/snakefmt/issues/115#issuecomment-986860293
if config["settings"].get("contig-group-size", 0) == 0:

    def contigs_groups_input(wildcards):
        return []


# Check that the config does not yet have contigs when this file is first included.
if "contigs" in config["global"]:
    raise Exception("Config key 'global:contigs' already defined. Someone messed with our setup.")


# Contigs in reference genome. This function gives the list of contig names that we want to
# use for calling, that is, either all contigs of the reference genome, or, if contig grouping
# is activated, the names for the contig groups that are sent as one job combining multiple contigs.
def get_contigs(fai):
    # This function might be called multiple times in different invocations of rules.
    # To avoid re-computing the contigs, we cache them in the global config dict.
    global config
    if "contigs" in config["global"]:
        return config["global"]["contigs"]

    # If the config sets a contig group size, we use this to solve a bin packing problem to
    # combine small contigs into a set, where each bin is at max as big as the threshold.
    # Here, we request the file via its checkpoit, to make sure that it is created by its rule
    # before we continue. This is valid, as this function here is only ever called from
    # within input functions of rules, which themselves request the fai file via checkpoint as well.
    if config["settings"].get("contig-group-size", 0) > 0:
        # Get the contigs group file. We parse it as a dict, whose keys are the contig group names.
        # Python wants us to explicitly convert this to a list here, as otherwise, some weird
        # pickling issue occurs downstream when snakemake tries to pickle the config for usage
        # in other rules...
        contig_group_file = checkpoints.contig_groups.get().output[0]
        contigs = json.load(open(contig_group_file))
        config["global"]["contigs"] = list(contigs.keys())
        return config["global"]["contigs"]

    # Without contig groups, just read the fai and return its first column,
    # which contains the ref sequence names (our contigs), except for the ones that are shorter
    # than the user-provided min size (if given). Store it in the global variable
    # first to not have to do the reading each time.
    min_contig_size = config["settings"].get("min-contig-size", 0)
    contig_list = read_contigs_from_fai(fai, min_contig_size)
    config["global"]["contigs"] = [t[0] for t in contig_list]
    return config["global"]["contigs"]


# =================================================================================================
#     Restrict Regions
# =================================================================================================

# Intersect the restict regions file with a given contig (chromosome), so that we can use the
# resulting bed file for parallelization over contigs.
if "restrict-regions" in config["settings"]:

    rule compose_regions:
        input:
            config["settings"].get("restrict-regions"),
        output:
            "calling/regions/{contig}.bed",
        log:
            "logs/calling/regions/bedextract/{contig}.regions.log",
        conda:
            "../envs/bedops.yaml"
        shell:
            "sort-bed {input} > {input}.sorted ; "
            "bedextract {wildcards.contig} {input}.sorted > {output}"

    # Rule is not submitted as a job to the cluster.
    localrules:
        compose_regions,


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
