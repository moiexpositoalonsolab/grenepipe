import json

# =================================================================================================
#     Get Fai
# =================================================================================================

def get_fai(wildcards):
    # Stop at the snakemake checkpoint first to ensure that the fai file is available.
    return checkpoints.samtools_faidx.get().output[0]
    # return config["data"]["reference"]["genome"] + ".fai"

# =================================================================================================
#     Grouping of (Small) Contigs
# =================================================================================================

if config["settings"].get("contig-group-size", 0) > 0:

    # Simple greedy solver for the bin packing problem, here used to make bins of roughly equal
    # size for the contigs. We except the input values to be tuples or pairs, where the index [1]
    # element is the weight (index[0] here is used for the contig name itself).
    def solve_bin_packing( values, max_bin_size ):
        # First, sort by length, decreasing.
        # This helps to get closer to an optimal solution.
        values.sort(key = lambda x: x[1], reverse=True)

        # Fill the bins as needed, using first-fit on the sorted list, and keeping track of
        # how much we already put in each of them. We can have at most as many bins as elements,
        # so use this to initialize the sums, to keep it simple.
        bins = []
        sums = [0] * len(values)
        for cont in values:
            # Find the first bin where the contig fits in.
            j = 0
            while( j < len(bins) ):
                if( sums[j] + cont[1] <= max_bin_size ):
                    bins[j].append( cont )
                    sums[j] += cont[1]
                    break
                j += 1

            # If no bin could fit the contig, make a new bin.
            # This in paticular catches contigs that are larger than the bin size, for example,
            # if some chromosomes are fully assembled and are hence not in scaffolds of small size.
            if j == len(bins):
                bins.append([])
                bins[j].append( cont )
                sums[j] += cont[1]

        # Log output
        # print("Contig group bin assignments for small contigs:")
        # for i in range(len(bins)):
        #     print(str(i) + ": [" + str(sums[i]) + "] " + str(bins[i]))

        # Return the bins with all contigs and their sizes
        return bins

    checkpoint contig_groups:
        input:
            fai = get_fai
        output:
            "contig-groups/contigs.json"
        log:
            "logs/contig-groups/contigs.log"
        params:
            contig_group_size = config["settings"].get("contig-group-size", 0)
        run:
            # Read fai to get all contigs and their sizes.
            contig_list = []
            with open(input.fai, "r") as f:
                for line in f:
                    contig, length_str = line.split("\t")[:2]
                    contig = contig.strip()
                    length = int(length_str.strip())
                    contig_list.append(( contig, length ))

            # Solve the bin packing for the contigs, to get a close to optimal solution
            # for putting them in groups. Large contigs (e.g., whole chromosomes) that are larger
            # than the bin size will simply get their own (overflowing...) bin.
            contig_bins = solve_bin_packing( contig_list, params.contig_group_size )

            # Now turn the contig bins into groups for the result of this function.
            # We store our resulting list of contigs containing all contigs,
            # in tuples with their sizes. The contigs is a dict from group name to a list
            # (over contings in the group) of tuples.
            contigs = {}
            for bin in contig_bins:
                groupname = "contig-group-" + str(len(contigs))
                contigs[groupname] = bin

            # We need to store the result in a file, so that the rule that creates the per-contig
            # files can access it.
            json.dump( contigs, open( output[0], 'w' ))

    # Rule is not submitted as a job to the cluster.
    localrules: contig_groups

    # Make the contig-group list files that contain the names of the contigs/scaffolds
    # that have been bin-packed above.
    rule contigs_group_list:
        input:
            contigs="contig-groups/contigs.json"
        output:
            "contig-groups/{contig}.bed"
        log:
            "logs/contig-groups/{contig}.log"
        run:
            # Get the contigs file that we created above.
            contigs = json.load( open( input.contigs ))

            # Same for the group name itself: This rule is only executed
            # for group names that we actually have made.
            if wildcards.contig not in contigs:
                raise Exception( "Internal error: contig " + wildcards.contig + " not found." )

            # Write the output list file, using the contig names and lengths from the group.
            # In bed files, the first three columns are the chrom name, start (incluse, zero-based),
            # and end (exclusive). Hence, we can simply use 0 and length as start and end here.
            with open(output[0], "w") as f:
                f.writelines( f"{c[0]}\t0\t{c[1]}\n" for c in contigs[wildcards.contig] )

    # Rule is not submitted as a job to the cluster.
    localrules: contigs_group_list

    # Conflicts of interest:
    if config["settings"].get("restrict-regions"):
        raise Exception(
            "Cannot combine settings contig-group-size > 0 with restrict-regions "
            "at the moment, as we have not implemented this yet. "
            "If you need this combination of settings, please submit an issue to "
            "https://github.com/lczech/grenepipe/issues and we will see what we can do."
        )

    # We now extended the rules for bcftools and freebayes to also work with small contig groups.
    # The following check is no longer needed - just kept here for reference.
    # if config["settings"]["calling-tool"] != "haplotypecaller":
    #     raise Exception(
    #         "Can only use setting contig-group-size with calling-tool haplotypecaller "
    #         "at the moment, as we have not implemented this for other calling tools yet. "
    #         "If you need this combination of settings, please submit an issue to "
    #         "https://github.com/lczech/grenepipe/issues and we will see what we can do."
    #     )

# =================================================================================================
#     Get Contigs
# =================================================================================================

# Check that the config does not yet have contigs when this file is first included.
if "contigs" in config["global"]:
    raise Exception("Config key 'global:contigs' already defined. Someone messed with our setup.")

# Contigs in reference genome.
def get_contigs( fai ):
    # If we have already computed contigs, just return them.
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
        contigs = json.load( open( contig_group_file ))
        config["global"]["contigs"] = list(contigs.keys())
        return config["global"]["contigs"]

    # Without contig groups, just read the fai and return its first column,
    # which contains the ref sequence names (our contigs). Store it in the global variable
    # first to not have to do the reading each time.
    config["global"]["contigs"] = pd.read_csv(
        fai, sep='\t', header=None, usecols=[0], squeeze=True, dtype=str
    )
    return config["global"]["contigs"]

# =================================================================================================
#     Restrict Regions
# =================================================================================================

# Intersect the restict regions file with a given contig (chromosome), so that we can use the
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
