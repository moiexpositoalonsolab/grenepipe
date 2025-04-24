Overview
============

Profiles that might come in handy as examples when running grenepipe locally or in a cluster setting. They are meant for the basic configuration, such as restart attempts, conda, etc. The profile in `slurm` also contains the basic slurm configuration of account and partition.

Note that the resource specifications for rule jobs are specified via the `config/resources.yaml` file since grenepipe v0.16.0, instead of specifying them here in the slurm config.

See the [Cluster and Profiles](https://github.com/lczech/grenepipe/wiki/Cluster-and-Profiles) wiki page for details on how those can be used with grenepipe. We also highly recommend to get familiar with the general Snakemake [Profiles])(https://snakemake.readthedocs.io/en/v8.15.2/executing/cli.html#profiles) as well as the Snakemake [SLURM Executor Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) if you want to run grenepipe on a cluster.
