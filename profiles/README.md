Overview
============

Profiles that might come in handy when running the pipeline.

The `local` profile is meant as a helper to not having to set all necessary settings by hand
each time snakemake is called.

The `slurm` profile is based on the snakemake profile cookiecutter template for slurm from
https://github.com/Snakemake-Profiles/slurm

However, we extended it as follows:
 - Slurm log files are collected in a subdirectory, instead of cluttering the main directory.
 - Using the `host` config files, specific configurations for each host can be provided.

This `host` config simply needs a file with the short hostname in the `host` subdirectory.
This file can contain any snakemake cluster configuration.
To get the short hostname on a given system, simply call `hostname -s`.

Usage
============

Direct usage example: `snakemake --profile profiles/slurm` for the `slurm` profile.

Alternatively, snakemake looks for profiles in `~/.config/snakemake`. Hence, you can also copy
the contents of the `local` or the `slurm` subdirectory to that location, and then do not need to
specify `--profile` when calling snakemake.
