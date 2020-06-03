Overview
============

Profiles that might come in handy when running the pipeline.

The `local` profile is meant as a helper to not having to set all necessary settings by hand
each time snakemake is called.

The `slurm` profile is based on the snakemake profile cookiecutter template for slurm from
https://github.com/Snakemake-Profiles/slurm

Usage example: `snakemake --profile path/to/local` for the `local` profile.
