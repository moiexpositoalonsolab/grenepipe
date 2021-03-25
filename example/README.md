Overview
==============

This directory contains some simple test/exemplary data, to get started with the pipeline.
See https://github.com/lczech/grenepipe/wiki/Quick-Start-and-Full-Example for details.


Samples
==============

Some data from the GrENE-net (http://grene-net.org) experiment, 10k random sequences in each sample.
Sample 2 has been sequenced twice, giving us two "units" to be analyzed in the pipeline.


Reference Genome
==============

Arabidopsis thaliana reference genome, obtained from

    https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas

The file `TAIR10_chr_all.fas.gz` has not to be decompressed, as this step is done automatically in the `prep` rules.


Pipeline
==============

To analyze this dataset, you first need to set the file paths in `config.yaml` and `samples.tsv` to fit with where
grenepipe is located, and then run

    snakemake --use-conda --cores 4 --directory example/ --snakefile rules/prep.smk
    snakemake --use-conda --cores 4 --directory example/

This calls variants and produces quality control statistics.
