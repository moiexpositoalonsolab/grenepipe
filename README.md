[![Snakemake](https://img.shields.io/badge/snakemake-8.15.2-brightgreen.svg)](https://snakemake.github.io)
[![CI](https://github.com/moiexpositoalonsolab/grenepipe/actions/workflows/ci.yaml/badge.svg)](https://github.com/moiexpositoalonsolab/grenepipe/actions)
[![Platforms](https://img.shields.io/badge/platform-linux--64%20%7C%20osx--64-lightgrey)](https://github.com/moiexpositoalonsolab/grenepipe/releases)
[![Release](https://img.shields.io/github/v/release/moiexpositoalonsolab/grenepipe.svg)](https://github.com/moiexpositoalonsolab/grenepipe/releases)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](http://www.gnu.org/licenses/gpl.html)
[![DOI](https://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbtac600-blue)](https://doi.org/10.1093/bioinformatics/btac600)
<!-- [![CI](https://github.com/moiexpositoalonsolab/grenepipe/workflows/CI/badge.svg?branch=master)](https://github.com/moiexpositoalonsolab/grenepipe/actions) -->
<!-- ![Language](https://img.shields.io/badge/language-python-lightgrey.svg) -->

![grenepipe logo](/doc/logo/grenepipe.png?raw=true)

Snakemake workflow for variant calling from raw sample sequences, with lots of bells and whistles.

**Advantages**:

  - One command to run the whole pipeline!
  - Many tools to choose from for each step
  - Simple configuration via a single file
  - Automatic download of tool dependencies
  - Resuming from failing jobs

Getting Started
-------------------

See [**--&gt; the Wiki pages &lt;--**](https://github.com/lczech/grenepipe/wiki) for setup and documentation.

For **questions, bug reports, and feature requests**,
[open an issue](https://github.com/lczech/grenepipe/issues).
Please do not send emails with questions or requests, as others might be having them as well,
and so it is better to discuss them where they can be found.

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
  - Read mapping
    - [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml)
    - [bwa aln](http://bio-bwa.sourceforge.net/bwa.shtml)
    - [bwa mem2](https://github.com/bwa-mem2/bwa-mem2)
    - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  - Optional read filtering, clipping, duplication removal, and quality score recalibration
    - [samtool view](http://www.htslib.org/doc/samtools-view.html)
    - [BamUtil clipOverlap](https://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap)
    - [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
    - [DeDup](https://github.com/apeltzer/dedup)
    - [GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator) (BQSR)
    - [samtools mpileup](http://www.htslib.org/doc/samtools-mpileup.html)
  - Damage profiling (optional; e.g., for ancient DNA)
    - [mapDamage](https://github.com/ginolhac/mapDamage)
    - [DamageProfiler](https://github.com/Integrative-Transcriptomics/DamageProfiler)
  - Variant calling and genotyping
    - [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) / [GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)
    - [freebayes](https://github.com/freebayes/freebayes)
    - [bcftools call](http://samtools.github.io/bcftools/bcftools.html#call)
  - Variant filtering
    - [GATK VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360036834871-VariantFiltration)
    - [GATK VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator) (VQSR)
    - [bcftools filter](https://samtools.github.io/bcftools/bcftools.html#filter)
  - Frequency calling (for pool sequencing data, as an alternative to variant calling)
    - [HAF-pipe](https://github.com/petrov-lab/HAFpipe-line)
  - Quality control, statistics, SNP annotation, reporting
    - [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    - [samtool stats](http://www.htslib.org/doc/samtools-stats.html)
    - [samtool flagstat](http://www.htslib.org/doc/samtools-flagstat.html)
    - [QualiMap](http://qualimap.conesalab.org/)
    - [Picard CollectMultipleMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360042478112-CollectMultipleMetrics-Picard-)
    - [bcftools stats](http://samtools.github.io/bcftools/bcftools.html#stats)
    - [snpEff](https://pcingola.github.io/SnpEff/)
    - [VEP](https://uswest.ensembl.org/info/docs/tools/vep/index.html) (Ensembl Variant Effect Predictor)
    - [MultiQC](https://multiqc.info/)

**Typical output:**

  - Variant calls `vcf`, raw and filtered, and potentially with annotations
  - MultiQC report (includes summaries of most other tools, and of the final `vcf`)
  - Snakemake report (optional)

Intermediate output files such as `bam` files are also kept by default,
and `mpileup` files can optionally be created if needed.
In addition to the above tools, there are some tools used as glue between the steps.
If you are interested in the details, have a look at the snakemake [rules](https://github.com/lczech/grenepipe/tree/master/rules) for each step.

Citation
-------------------

When using grenepipe, please cite:

> **grenepipe: A flexible, scalable, and reproducible pipeline <br/>to automate variant calling from sequence reads.**<br/>
> Lucas Czech and Moises Exposito-Alonso. *Bioinformatics*. 2022.<br/>
> [doi:10.1093/bioinformatics/btac600](https://doi.org/10.1093/bioinformatics/btac600) [[pdf](https://drive.google.com/file/d/125IRw_orGGxWWYr5GZ1LMCHbFDXj0C04/view?usp=sharing)]

Furthermore, please do not forget to cite all tools that you selected to be run for your analysis. See [our Wiki](https://github.com/moiexpositoalonsolab/grenepipe/wiki/Citation-and-References) for their references.
