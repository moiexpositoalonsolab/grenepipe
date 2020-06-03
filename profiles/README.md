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

Usage example: `snakemake --profile profiles/slurm` for the `slurm` profile.
