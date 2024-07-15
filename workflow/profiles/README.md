Overview
============

Profiles that might come in handy when running the pipeline in a cluster setting. The profile in `slurm` also contains a basic slurm configuration for some of the rule time and memory requirements that have worked for us for variant calling on normal-sized fastq inputs.

See the [Cluster and Profiles](https://github.com/lczech/grenepipe/wiki/Cluster-and-Profiles) wiki page for details on how those can be used with grenepipe. We also highly recommend to get familiar with the general Snakemake [Profiles])(https://snakemake.readthedocs.io/en/v8.15.2/executing/cli.html#profiles) as well as the Snakemake [SLURM Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) if you want to run grenepipe on a cluster.
