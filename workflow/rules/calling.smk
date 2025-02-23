import json
import math

# =================================================================================================
#     Get Contigs
# =================================================================================================


# Get the list of chromosome names that are present in the fai file,
# and their length, with a length filter if needed.
def read_contigs_from_fai(fai, min_contig_size = 0):
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


# =================================================================================================
#     Grouping of (Small) Contigs
# =================================================================================================

if config["settings"].get("contig-group-size", 0) > 0:

    # Simple greedy solver for the bin packing problem, here used to make bins of roughly equal
    # size for the contigs. We expect the input values to be tuples or pairs, where the index [1]
    # element is the weight that we use for packing (index[0] here is used for the contig name itself).
    def solve_bin_packing(values, max_bin_size):
        # First, sort by length, large ones first, decreasing.
        # This helps to get closer to an optimal solution.
        values.sort(key=lambda x: x[1], reverse=True)

        # Fill the bins as needed, using first-fit on the sorted list, and keeping track of
        # how much we already put in each of them. We can have at most as many bins as elements,
        # so use this to initialize the sums, to keep it simple.
        bins = []
        sums = [0] * len(values)
        for cont in values:
            # Find the first bin where the contig fits in.
            j = 0
            while j < len(bins):
                if sums[j] + cont[1] <= max_bin_size:
                    bins[j].append(cont)
                    sums[j] += cont[1]
                    break
                j += 1

            # If no bin could fit the contig, make a new bin.
            # This in paticular catches contigs that are larger than the bin size, for example,
            # if some chromosomes are fully assembled and are hence not in scaffolds of small size.
            if j == len(bins):
                bins.append([])
                bins[j].append(cont)
                sums[j] += cont[1]

        # Log output
        # print("Contig group bin assignments for small contigs:")
        # for i in range(len(bins)):
        #     print(str(i) + ": [" + str(sums[i]) + "] " + str(bins[i]))

        # Return the bins with all contigs and their sizes
        print("Success with", len(bins), "groups")
        return bins


    # Alternative heuristic to assign contigs to groups, minimizing the number of groups created,
    # while also respecting the max contig group size and the max contigs per group, and also
    # trying to optimize the load balancing by distributing short and long contigs over groups.
    # The result is a list of lists of tuples. The outer list are the groups. Each group then
    # contains the list of tuples (contig name and length) of the group.
    def optimize_contig_group(contigs, max_contig_group_size, max_contigs_per_group):
        # First, sort by length, large ones first, decreasing.
        contigs.sort(key=lambda x: x[1], reverse=True)

        # First, fill the final group list with all contigs that are larger than the group size.
        # We do not split existing chromsomes/contigs, so we have to live with them being too large.
        # Also, collect all smaller contigs (smaller than the max group size), for further processing.
        large_groups = []
        small_contigs = []
        for contig in contigs:
            if contig[1] >= max_contig_group_size:
                large_groups.append([])
                large_groups[-1].append(contig)
            else:
                small_contigs.append(contig)
        print("Found ", len(large_groups), "large contigs and", len(small_contigs), "small contigs")

        # No small contigs: We are done
        if len(small_contigs) == 0:
            return large_groups

        # Now we have a list of small contigs for which to create new groups.
        # We do that in a loop until we succeed to fulfill all boundary conditions.
        # In particular, we want to limit the number of contigs per group, and not exceed
        # the sum of contig lengths per group. In each iteration, we increase the number
        # of groups created by one, hence guaranteeing success at some point.
        # There might be a better solution to this, but this is simple and works.
        # We start with as many groups as we minimally need such that we can fulfill the max
        # contigs per group constraint.
        num_groups = max(1, int(math.ceil(len(small_contigs) / max_contigs_per_group)))
        need_iteration = True
        while need_iteration:
            print("Evaluating with", num_groups, " small contig groups")
            need_iteration = False

            # We fill a temporary list for the small contigs, starting with as many empty lists
            # as we have groups, and then fill each group's list with its contigs.
            # We do this by alternating to fill forwards and backwards - that is the heuristic
            # that is meant to avoid having groups with only small or only large contigs.
            # By going back and forth instead of filling from the start, we get a more even spread.
            # Let's call this the zig-zag round-robing assignment :-)
            # For instance, we should get: [[1, 6, 7], [2, 5, 8], [3, 4, 9]]
            small_groups = [[] for _ in range(num_groups)]
            round_num = 0
            index = 0

            # Continue until all contigs have been distributed.
            while index < len(small_contigs):
                # Determine the order based on whether the round is even (forward) or odd (backward)
                if round_num % 2 == 0:
                    order = range(num_groups)
                else:
                    order = range(num_groups - 1, -1, -1)

                # Distribute one element per group in the specified order.
                for i in order:
                    if index >= len(small_contigs):
                        break
                    small_groups[i].append(small_contigs[index])
                    index += 1
                round_num += 1

            # Now we test if this was successful: Is each group small enough, both in terms
            # of total length, and in terms of number of contigs per group?
            small_contig_count = 0
            for group in small_groups:
                total_len = sum(contig[1] for contig in group)
                if (total_len > max_contig_group_size) or (len(group) > max_contigs_per_group):
                    print(
                        "Not enough groups. Got group length", total_len, ">", max_contig_group_size,
                        "and got contigs per group", len(group), ">", max_contigs_per_group
                    )
                    need_iteration = True
                    num_groups += 1
                    break
                small_contig_count += len(group)

        # Now we are done. Make sure that we have processed the right number of contigs.
        # Then, add all small contigs as groups to our final result, and return it.
        assert small_contig_count == len(small_contigs)
        print(
            "Success with", num_groups, "small contig groups, and",
            (len(large_groups) + len(small_groups)), "total groups"
        )
        return large_groups + small_groups


    checkpoint contig_groups:
        input:
            fai=get_fai,
        output:
            "calling/contig-groups/contigs.json",
        log:
            "logs/calling/contig-groups/contigs.log",
        params:
            min_contig_size=config["settings"].get("min-contig-size", 0),
            contig_group_size=config["settings"].get("contig-group-size", 0),
            max_contigs_per_group=config["settings"].get("max-contigs-per-group", 0),
        run:
            # Open the log file in write (or append) mode, and redirect both stdout and stderr
            # to the log file. Snakemake does not do that for us...
            log_file = open(log[0], "w")
            sys.stdout = log_file
            sys.stderr = log_file

            # Solve the bin packing for the contigs, to get a close to optimal solution
            # for putting them in groups. Large contigs (e.g., whole chromosomes) that are larger
            # than the bin size will simply get their own (overflowing...) bin.
            contig_list = read_contigs_from_fai(input.fai, params.min_contig_size)
            if params.max_contigs_per_group == 0:
                print("Running bin packing solver")
                contig_groups = solve_bin_packing(contig_list, params.contig_group_size)
            else:
                print("Running heuristic optimizer")
                contig_groups = optimize_contig_group(
                    contig_list, params.contig_group_size, params.max_contigs_per_group
                )

            # Now turn the contig bins into groups for the result of this function.
            # We store our resulting list of contigs containing all contigs,
            # in tuples with their sizes. The contigs is a dict from group name to a list
            # (over contings in the group) of tuples.
            contigs = {}
            for group in contig_groups:
                groupname = "contig-group-" + str(len(contigs))
                contigs[groupname] = group

            # We need to store the result in a file, so that the rule that creates the per-contig
            # files can access it.
            json.dump(contigs, open(output[0], "w"))


    # Rule is not submitted as a job to the cluster.
    localrules:
        contig_groups,

    # Due to a new bug in Snakemake after our update to 8.15.2, we now need the following
    # function to be called as input whenever the above checkpoint is needed.
    # See https://github.com/snakemake/snakemake/issues/2958 for details.
    def contigs_groups_input(wildcards):
        return checkpoints.contig_groups.get().output[0]

        # Make the contig-group list files that contain the names of the contigs/scaffolds
        # that have been bin-packed above.
    rule contigs_group_list:
        input:
            contigs="calling/contig-groups/contigs.json",
        output:
            "calling/contig-groups/{contig}.bed",
        # log:
        #     "logs/calling/contig-groups/{contig}.log",
        run:
            # Get the contigs file that we created above.
            contigs = json.load(open(input.contigs))

            # Same for the group name itself: This rule is only executed
            # for group names that we actually have made.
            if wildcards.contig not in contigs:
                raise Exception("Internal error: contig " + wildcards.contig + " not found.")

                # Write the output list file, using the contig names and lengths from the group.
                # In bed files, the first three columns are the chrom name, start (incluse, zero-based),
                # and end (exclusive). Hence, we can simply use 0 and length as start and end here.
            with open(output[0], "w") as f:
                f.writelines(f"{c[0]}\t0\t{c[1]}\n" for c in contigs[wildcards.contig])

                # Rule is not submitted as a job to the cluster.


    localrules:
        contigs_group_list,

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


# Dummy definition of the above rule for when we are not using contig groups.
# We cannot do this as an else branch of the above, due to an issue in the
# formatting of snakemfmt: https://github.com/snakemake/snakefmt/issues/115#issuecomment-986860293
if config["settings"].get("contig-group-size", 0) == 0:

    def contigs_groups_input(wildcards):
        return []


# =================================================================================================
#     Get Contigs
# =================================================================================================

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
    [t[0] for t in tuples_list]
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
