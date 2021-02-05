# =================================================================================================
#     README
# =================================================================================================

# Snakemake wrapper for freebayes, adapted from the wrapper script at
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/freebayes.html in version 0.70.0,
# https://github.com/snakemake/snakemake-wrappers/blob/52c4a4a4a8fb388e68b46c0a9982432ef86cd071/bio/freebayes/wrapper.py
#
# We use this as a basis for developing our own improved freebayes script
# that fixes a bug where with 1 thread, the regions would not be taken into account,
# and also has better scalability and load balancing / distribution on the cluster.
# --> As of now, we did not entirely succeed in the latter, as using `bamtools coverage` for getting
# equal workloads (a suggestion of freebayes itself) does not work properly... Bionformatics tools...
# But still, we get per-contig parallelization, which is something!

# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell

shell.executable("bash")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

params = snakemake.params.get("extra", "")
norm = snakemake.params.get("normalize", False)
assert norm in [True, False]

# Additions to the original wrapper made by LC:
#
# We separate by contig, to increase parallel runs on clusters.
# That is, the original wrapper just calls one instance of freebayes for the whole input,
# and runs parallel threads on the same compute node (in cluster environments) for different regions.
#
# This however still means that we are only using one compute node for the whole genome.
# We do better here, and instead split by contig, so that at least that many nodes can be active
# at the same time. Within the contigs, we then again split into regions, so that each node can
# also use multiple threads. But different from the original wrapper, our region split is not simpy
# using a fixed region length, but instead looks for actual content to ensure better parallel
# efficiency (threads that have a similar runtime).
#
# In theory, it would be better to have more fine grained ways of splitting, so that we could split
# the load to n compute nodes, instead of using contigs to split for the nodes.
# Future work...
contig = snakemake.wildcards.contig

# Get the fai for the ref genome. This is expected to follow a file naming convention.
# We could hence instead get this file from just the ref file name as well, as done in the
# original wrapper script, but by instead requiring it as a rule input, we make sure that it
# actually exists and get a nice snakemake error message if not.
fai = snakemake.input.fai
assert snakemake.input.ref + ".fai" == fai

# =================================================================================================
#     Main Call
# =================================================================================================

# Prepare the output and compression.
pipe = ""
if snakemake.output[0].endswith(".bcf"):
    if norm:
        pipe = "| bcftools norm -Ob -"
    else:
        pipe = "| bcftools view -Ob -"
elif norm:
    pipe = "| bcftools norm -"

# Get the full chromosome as a region in a bed compatible format for the given contig.
regions = ""
with open(fai) as faif:
    for line in faif:
        fields = line.strip().split("\t")
        chrom_name = fields[0]
        chrom_length = int(fields[1])
        if chrom_name == contig:
            regions = "<(echo \"" + chrom_name + ":0-" + str(chrom_length) + "\")"

# If we are here, we must have found the contig in the fai file,
# otherwise that name would not have appeared in the "{contig}" wildcard of our snakemake rule.
assert regions != ""

# Now intersect with the regions file (if provided). This is the bug fix compared to the original.
if snakemake.input.get("regions", ""):
    regions = (
        "<(bedtools intersect -a "
        r"<(sed 's/:\([0-9]*\)-\([0-9]*\)$/\t\1\t\2/' "
        "{regions}) -b {snakemake.input.regions} | "
        r"sed 's/\t\([0-9]*\)\t\([0-9]*\)$/:\1-\2/')"
    ).format(regions=regions, snakemake=snakemake)

if snakemake.threads == 1:
    freebayes = "freebayes --region <("+ regions + ")"
else:
    # Ideally, we'd be using bamtools coverage and coverage_to_regions.py here,
    # as suggsted in the freebayes-parallel script, but this runs a long time and had some errors
    # in our test data already, so for now, we just split into chunks of equal size again,
    # similar to what the original wrapper does.
    chunksize = snakemake.params.get("chunksize", 100000)
    chunks = "<(fasta_generate_regions.py {snakemake.input.ref}.fai {chunksize})".format(
        snakemake=snakemake, chunksize=chunksize
    )

    # Now we need to intersect those with the regions that we got from above.
    regions = (
        "<(bedtools intersect -a "
        r"<(sed 's/:\([0-9]*\)-\([0-9]*\)$/\t\1\t\2/' "
        "{regions}) -b "
        r"<(sed 's/:\([0-9]*\)-\([0-9]*\)$/\t\1\t\2/' "
        "{chunks}) | "
        r"sed 's/\t\([0-9]*\)\t\([0-9]*\)$/:\1-\2/')"
    ).format(regions=regions, chunks=chunks)
    freebayes = ("freebayes-parallel {regions} {snakemake.threads}").format(
        snakemake=snakemake, regions=regions
    )

    # Somehow, in our script here the call to `vcfstreamsort` that is done in the freebayes-parallel
    # script does not work. In fact, I don't understand how it can work for the original wrapper
    # script at all: It uses a sorting window with that is below the region width that is used
    # by default in the fasta regions script. So that should never fully sort everything.
    #
    # We could start modifying that script as well and provide our own version with a wider sorting
    # window width... But instead, let's just sort again here, with a better window size.
    #
    # If this still does not work at some point (because the freebayes calles are still too mixed),
    # we could try to increase the size to chunksize*threads (I think that is the maximum that
    # can happen in terms of unorder), or simply replace the following by a call to `| bcftools sort - `
    pipe = "| vcfstreamsort -w {chunksize} | vcfuniq ".format(chunksize=chunksize) + pipe

shell(
    "({freebayes} {params} -f {snakemake.input.ref}"
    " {snakemake.input.samples} {pipe} > {snakemake.output[0]}) {log}"
)
