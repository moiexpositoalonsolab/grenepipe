Overview
==============

This directory contains some simple test data, to get started with the pipeline.

Samples
==============

Some sampling data from garden experiment, 10k random sequences in each sample.
We artificially use sample B twice, once as a single end read, once as pair end reads,
just to make sure that the pipeline works with both cases.

Reference Genome
==============

Arabidopsis thaliana reference genome, obtained from

    https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

The file `genome.fa.gz` has not to be decompressed, as this step is done automatically in the `prep` rules.
