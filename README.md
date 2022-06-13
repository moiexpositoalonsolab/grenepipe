![grenepipe logo](/doc/logo/grenepipe.png?raw=true)

Snakemake pipeline for variant calling from raw sample sequences, with lots of bells and whistles.

**Advantages**:

  - One command to run the whole pipeline!
  - Many tools to choose from for each step
  - Simple configuration via a single file
  - Automatic download of tool dependencies
  - Resuming from failing jobs

Pipeline Overview
-------------------

**Minimal input:**

  - Reference genome `fasta` file
  - Per-sample `fastq` files
  - Optionally, a `vcf` file of known variants to restrict the variant calling process

**Process and available tools:**

  - Read trimming (single or paired end)
    - [AdapterRemoval](https://adapterremoval.readthedocs.io/en/latest/)
    - [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
    - [fastp](https://github.com/OpenGene/fastp)
    - [SeqPrep](https://github.com/jstjohn/SeqPrep)
    - [skewer](https://github.com/relipmoc/skewer)
    - [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
  - Read mapping, optionally with duplication removal and quality score recalibration
    - [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml)
    - [bwa aln](http://bio-bwa.sourceforge.net/bwa.shtml)
    - [bwa mem2](https://github.com/bwa-mem2/bwa-mem2)
    - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
    - [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
    - [DeDup](https://github.com/apeltzer/dedup)
    - [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator)
    - [samtools mpileup](http://www.htslib.org/doc/samtools-mpileup.html) (optional)
  - Damage profiling (optional; e.g., for ancient DNA)
    - [mapDamage](https://github.com/ginolhac/mapDamage)
    - [DamageProfiler](https://github.com/Integrative-Transcriptomics/DamageProfiler)
  - Variant calling, genotyping, and filtering
    - [bcftools call](http://samtools.github.io/bcftools/bcftools.html)
    - [freebayes](https://github.com/freebayes/freebayes)
    - [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) / [GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)
    - [GATK SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants)
    - [GATK VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360036834871-VariantFiltration)
    - [GATK VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator)
  - Quality control, statistics, SNP annotation, reporting
    - [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    - [samtool stats](http://www.htslib.org/doc/samtools-stats.html)
    - [samtool flagstat](http://www.htslib.org/doc/samtools-flagstat.html)
    - [QualiMap](http://qualimap.conesalab.org/)
    - [Picard CollectMultipleMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360042478112-CollectMultipleMetrics-Picard-)
    - [snpEff](https://pcingola.github.io/SnpEff/)
    - [VEP](https://uswest.ensembl.org/info/docs/tools/vep/index.html) (Ensembl Variant Effect Predictor)
    - [MultiQC](https://multiqc.info/)

**Output:**

  - Variant calls `vcf`, and table of allele frequencies (optional, e.g., for pool sequencing data)
  - MultiQC report (includes summaries of most other tools)
  - Snakemake report (optional)

Intermediate output files such as `bam` files are also kept,
and `mpileup` files can optionally be created if needed.
In addition to the above tools, there are some tools used as glue between the steps.
If you are interested in the details, have a look at the snakemake [rules](https://github.com/lczech/grenepipe/tree/master/rules) for each step.

Getting Started
-------------------

See [**the Wiki pages**](https://github.com/lczech/grenepipe/wiki) for setup and documentation.

For **bug reports and feature requests**, please
[open an issue on our GitHub page](https://github.com/lczech/grenepipe/issues).

Citation
-------------------

When using grenepipe, please cite:

> **grenepipe: A flexible, scalable, and reproducible pipeline to automate variant and frequency calling from sequence reads.**<br/>
> Lucas Czech and Moises Exposito-Alonso. *arXiv*. 2021.<br/>
> [doi:10.48550/arXiv.2103.15167](https://doi.org/10.48550/arXiv.2103.15167), [arXiv:2103.15167](https://arxiv.org/abs/2103.15167)

Furthermore, please do not forget to cite all tools that you selected to be run for your analysis. See [our Wiki](https://github.com/moiexpositoalonsolab/grenepipe/wiki/Citation-and-References) for their references.
