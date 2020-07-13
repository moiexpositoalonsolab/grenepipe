# grenepipe

Snakemake pipeline for allele frequency calculations in evolve &amp; resequence experiments.
The pipeline roughly follows the [GATK best practices workflow](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows), but with a lot of additional options and tool choices.

Input:
  - Reference genome
  - Per-sample `fastq` files, either single end or paired end, and potentially with multiple units,
    that is, different sequencing runs or multiple lanes per sample
  - Optionally, a vcf file of known variants to restrict the variant calling process.

Process:
  - Read trimming
  - Read mapping, optionally with duplication removal and quality score recalibration
  - Variant calling, genotyping, and filtering
  - Quality control, SNP annotation

Output:
  - Variant calls `vcf`

## Basic Usage

The file `config.yaml` contains the basic configuration for the data and the tools to use:

 *  Path to a list of all sample `fastq` files to process, see below for the expected format.
 *  Path to the reference genome `fasta` file.
 *  Settings for which tools to use and their parameters.

Change this file to your needs.

### Samples Table

The configuration (`config.yaml`) expects the key `data --> samples` to point to a tab-separated
table that lists all sample names and the paths to their `fastq` files:

    sample	unit	platform	fq1	fq2
    A	1	ILLUMINA	data/samples/A_R1_001.fastq.gz	data/samples/A_R2_001.fastq.gz
    B	1	ILLUMINA	data/samples/B_R1_001.fastq.gz	data/samples/B_R2_001.fastq.gz
    B	2	ILLUMINA	data/samples/B_002.fastq.gz

This table contains two samples, `A` and `B`. Sample `A` consists of one unit (coming from one
sequencing run), with paired-end reads. Sample `B` consists of two units (coming e.g., from
different sequencing runs or multiple lanes per sample), with the first unit being paired-end,
and the second unit being single-end (hence, the last column is simply left empty).
Samples can either be in `fastq` format, or in compressed `fastq.gz` format (the file extension
itself does not matter though).

### Reference Genome

The reference genome is expected to be in `fasta` format. Furthermore, we do annotation using
[snpEff](http://snpeff.sourceforge.net/), which expects a valid genome name from its database
in our `config.yaml` under the key `data --> reference --> name`; see below for details.

### Preparing the Pipeline

There are some steps that need to be done before running the main data analysis.
These are mainly for preparing indices into the reference genome, which decide how our jobs are
split up. In order for Snakemake to use this information, these steps (currently) have to be run
before starting the actual analysis:

    snakemake --cores 6 --use-conda --snakefile rules/prep.smk

will run the `prep` rules of the pipeline.

### Running the Pipeline

Once this is done, simply run

    snakemake --cores 6 --use-conda

for the whole pipeline.

### Advanced Usage

The above call will produce output files in the main directory. This might get crowded.
Also, one might want to use different configurations, i.e., different sets of input data,
different choices of tools and their parameters, etc.

Hence, the snakemake `--directory` can be set, from which the config file is read and to which
output files are written:

    --directory my-test

Caveat: With that setup, snakemake unfortunately stores all conda environments inside the specified
directory. That means, all conda environments are downloaded again and again for each `--directory`
that we use. This is of course not desirable (and it is mysterious why snakemake behaves that way),
so let's avoid that by setting

    --conda-prefix /path/to/some/stable/conda/directory

This setup also allows to run snakemake with different datasets at the same time, e.g,
in a cluster environment.

### Cluster Environments

To fully leverage cluster environments, rules can be executed as individual jobs in parallel.
Snakemake provides so-called [profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
to wrap the configurations necessary for each cluster enviroment.

For our purposes, we have adapted and refined the [slurm profile](https://github.com/Snakemake-Profiles/slurm),
which can be activated by calling

    --profile profiles/slurm

See the [profiles readme](profiles/README.md) for details.

### Report

We also offer to automatically generate a Snakemake report for a run of the pipeline

    snakemake --directory my-test --report my-test/report.html

This can of course also be used without the `--directory` option.

For this to work, the Python packages networkx and pygraphviz must be installed:

    sudo apt-get install python3-dev graphviz libgraphviz-dev pkg-config python-pip
    sudo pip install networkx pygraphviz

## Rule Dependencies and Call Graph

The default setup, using `bwa mem` for mapping, and `GATK HaplotypeCaller` for calling,
has the following rule dependencies:

![Rule dependency graph](/tools/rules.png)

If configured to use different tools, such as `bowtie2` and `bcftools call`,
the graph might look different, depending on which extra steps are needed for these tools.

To get the rule dependency graph and the job call graph, use

    snakemake --rulegraph | dot -Tsvg > rules.svg
    snakemake --dag | dot -Tsvg > dag.svg

## snpEff

The annotation tool snpEff needs a bit of special consideration here.
For it to work, we need to select a reference genome that snpEff understands.

1.  First, we install snpEff in a conda environment, so that we can work cleanly.

    Create the file `env-snpeff.yaml`:
    ```
    channels:
      - conda-forge
      - bioconda
    dependencies:
      - snpeff =4.3.1t
    ```

    Now, create and activate the conda environment:
    ```
    conda env create --name snpeff --file env-snpeff.yaml
    conda activate snpeff
    ```

2.  Next, find the name of your reference genome in the snpEff database.
    Here, we filter the database for some term that we expect to find.
    Of course, change "thaliana" to what you are looking for,
    or just omit the grep to get the whole list.

    ```
    snpEff databases | grep -i "thaliana"
    ```

2.  The first column of the output is what we are looking for:
    `Arabidopsis_thaliana` is the name of the reference genome for A. thaliana in the snpEff database.
    Put this name into the `config.yaml` file at the `data --> reference --> name` key.
