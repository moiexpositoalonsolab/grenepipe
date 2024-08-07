Overview
==============

This directory contains some simple exemplary and testing data, to get started with the pipeline.
See https://github.com/lczech/grenepipe/wiki/Quick-Start-and-Full-Example for details.


Samples
==============

Some data from the GrENE-net (http://grene-net.org) experiment, 10k random sequences in each sample.
Sample 2 has been sequenced twice, giving us two "units" to be analyzed in the pipeline.
For sample 3, we only used the R1 strand, giving us a (fake) example of single-end file for testing,
although for completeness, we include both the R1 and R2 files here (R2 is just not listed
in the samples table).

The original source files from our data are as follows:

 * Sample 1: `MLFH10_1_20180130` (data release 1)
 * Sample 2: `MLFH_9_1_20180531` (data release 1 and 3)
 * Sample 3: `MLFH_9_2_20180531` (data release 3)


Reference Genome
==============

Arabidopsis thaliana reference genome, obtained from

    https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

The file `TAIR10_chr_all.fas.gz` has not to be decompressed, as this step is done automatically.


Known Variants and Regions
==============

For testing purposes, we here also include a `known-variants.vcf.gz` file that can be used
to constrain the variant calling process. This file is based on the [1001 Genomes](https://www.1001genomes.org/) dataset,
and subset to serve for exemplary and test purposes (see `known-variants.sh` for the subsetting steps if you are interested).

Furthermore, for testing, we provide a simple `regions.bed` file that restricts the calling process
to just some of the regions in the chromosomes. Note that this is mostly an experimental feature
at the moment (see the `config.yaml`, under `settings: restrict-regions`, for details).

Lastly, the `mapping` directory contains mapped `bam` files of the samples, also for testing
purposes. This can for instance be used in order to run the pipeline for variant calling
starting from a set of given `bam` files, instead of mapping `fastq` files first.
See [here](https://github.com/moiexpositoalonsolab/grenepipe/wiki/Advanced-Usage#running-only-parts-of-the-pipeline) for details on this.


Pipeline
==============

To analyze this exemplary dataset, follow these steps:

 * First, copy the file `./config/config.yaml` into this example directory.
 * Then set the file paths in `config.yaml` and `samples.tsv` to match where this `example` directory is located.

Alternative to these two steps, simply call `prepare.sh`, which does these steps automatically.

Finally, from the main grenepipe directory, run

    snakemake --use-conda --conda-frontend mamba --cores 4 --directory example/

This calls variants and produces quality control statistics.
See https://github.com/lczech/grenepipe/wiki/Quick-Start-and-Full-Example for details.
