import json

# =================================================================================================
#     Common Helper Functions
# =================================================================================================

def get_fai(wildcards):
    # Stop at the snakemake checkpoint first to ensure that the fai file is available.
    return checkpoints.samtools_faidx.get().output[0]
    # return config["data"]["reference"]["genome"] + ".fai"

# Check that the config does not yet have contigs when this file is first included.
if "contigs" in config["global"]:
    raise Exception("Config key 'global:contigs' already defined. Someone messed with our setup.")

def solve_contig_bin_packing( small_contigs, small_contig_thresh ):
    # Sort by longest (of the small) contig first.
    # This helps to get closer to an optimal solution.
    small_contigs.sort(key = lambda x: x[1], reverse=True)

    # Fill the bins as needed, using first-fit on the sorted list, and keeping track of
    # how much we already put in each of them. We can have at most as many bins as contigs.
    bins = []
    sums = [0] * len(small_contigs)
    for cont in small_contigs:
        # Find the first bin where the contig fits in.
        j = 0
        while( j < len(bins) ):
            if( sums[j] + cont[1] <= small_contig_thresh ):
                bins[j].append( cont )
                sums[j] += cont[1]
                break
            j += 1

        # If no bin could fit the contig, make a new bin.
        if j == len(bins):
            bins.append([])
            bins[j].append( cont )
            sums[j] += cont[1]

    # Debug print
    # for i in range(len(bins)):
        # logger.info(str(i) + ": [" + str(sums[i]) + "] " + str(bins[i]))

    # Return the bins with all contigs and their sizes
    return bins

def make_small_contig_groups( fai ):
    global config
    small_contig_thresh = config["settings"].get("small-contigs-threshold", 0)
    assert small_contig_thresh > 0

    # We store our resulting list of contig names, as well as a list of bins,
    # containing all (large and small) contigs, in tuples with their sizes.
    # The contig-bins is a dict from group name to a list (over contings in the group) of tuples.
    assert "contigs" not in config["global"]
    config["global"]["contigs"] = []
    config["global"]["contig-bins"] = {}
    contig_cnt = 0
    small_contigs = []

    # Read fai to get all contigs and their sizes.
    # Put the large ones into the result immediately, and collect the small ones in a list
    # of pairs, with name and size of each contig, so that we can run the bin packing.
    with open(fai, "r") as f:
        for line in f:
            contig, length_str = line.split("\t")[:2]
            contig = contig.strip()
            length = int(length_str.strip())

            if length >= small_contig_thresh:
                # Large ones are immediately added to the result.
                groupname = "contig-group-" + str(contig_cnt)
                config["global"]["contigs"].append(groupname)
                config["global"]["contig-bins"][groupname] = [( contig, length )]
                contig_cnt += 1
            else:
                small_contigs.append(( contig, length ))

    # Solve the bin packing for the small contigs, to get a close to optimal solution
    # for putting them in groups.
    bins = solve_contig_bin_packing( small_contigs, small_contig_thresh )

    # Now turn the small contig bins into groups for the result of this function.
    for bin in bins:
        groupname = "contig-group-" + str(contig_cnt)
        config["global"]["contigs"].append(groupname)
        config["global"]["contig-bins"][groupname] = bin
        contig_cnt += 1

    # Debug
    # logger.info(str(config["global"]["contig-bins"]))

    # We need to store the result in a file, so that the rule that creates the per-contig files
    # can access it. This is super hacky, as we are currently not in a rule here, so technically,
    # we should not rely on file paths etc here... But it seems that snakemake sets the working
    # directory globally, so that our file access here puts the data in the right place.
    # Alternative solutions might be a persistent dict, but that requires additionaly python
    # modules to be imported, or to call this function again from a rule, but that would
    # require to solve the bin packing twice. The algorithm is deterministic, so that should not
    # change the result, but still it seems not like a good solution either, so let's roll with
    # this one here and hope that the working directory is always correct.
    os.makedirs("contig-groups", exist_ok=True)
    json.dump( config["global"]["contig-bins"], open( "contig-groups/contigs.json", 'w' ))

# Contigs in reference genome.
def get_contigs( fai ):
    # If we have already computed contigs, just return them.
    global config
    if "contigs" in config["global"]:
        return config["global"]["contigs"]

    # If the config sets a small contig threshold, we use this
    # to solve a bin packing problem to combine small contigs into a set, where each bin
    # is at max as big as the threshold.
    small_contig_thresh = config["settings"].get("small-contigs-threshold", 0)
    if small_contig_thresh > 0:
        make_small_contig_groups( fai )
        assert "contigs" in config["global"]
        return config["global"]["contigs"]

    # Without small contig threshold, just read the fai and return its first column,
    # which contains the ref sequence names (our contigs). Store it in the global variable
    # first to not have to do the reading each time.
    config["global"]["contigs"] = pd.read_csv(
        fai, sep='\t', header=None, usecols=[0], squeeze=True, dtype=str
    )
    return config["global"]["contigs"]

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

    # Make the contig-group list files that contain the names of the contigs/scaffolds
    # that have been bin-packed above.
    rule group_contigs:
        input:
            ref=get_fai
        output:
            "contig-groups/{contig}.list"
        log:
            "logs/contig-groups/{contig}.log"
        run:
            # Get the contigs file that we created above. This will always exist, as this rule
            # is only ever executed after we have resolved the "{contig}" var, which means,
            # get_contigs() was already called.
            if not os.path.exists( "contig-groups/contigs.json" ):
                raise Exception( "Internal error: get_contigs() was not called yet." )
            contigs = json.load( open( "contig-groups/contigs.json" ))

            # Same for the group name itself: This rule is only executed
            # for group names that we actually have made.
            if wildcards.contig not in contigs:
                raise Exception( "Internal error: contig " + wildcards.contig + " not found." )

            # Write the output list file, using the contig names from the group.
            with open(output[0], "w") as f:
                f.writelines( f"{c[0]}\n" for c in contigs[wildcards.contig] )

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
