# Configuring grenepipe

Grenepipe is a highly flexible workflow for variant calling from raw sample sequences,
with lots of bells and whistles. To configure this workflow, modify `config/config.yaml`
according to your needs, following the explanations provided in the file.

Furthermore, for the general usage of grenepipe, see our
[wiki](https://github.com/moiexpositoalonsolab/grenepipe/wiki).
See there to get started with grenepipe.

# Pipeline Overview

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
