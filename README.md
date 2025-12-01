# RNA-seq Processing Pipeline (Nextflow DSL2)

This repository contains a comprehensive **container-less RNA-seq
analysis pipeline** built using **Nextflow DSL2**.\
It supports raw FASTQ input and performs quality control, trimming,
alignment, quantification, UMI processing, and various QC/metrics steps
using widely used bioinformatics tools.

------------------------------------------------------------------------

## 📌 Features

The pipeline includes:

### Quality Control

-   **FastQC** --- initial raw read QC\
-   **MultiQC** --- aggregated QC report

### Pre-processing

-   **Fastp** --- adapter trimming & filtering\
-   **UMI-tools** --- UMI extraction and deduplication\
-   **SortMeRNA** (optional) --- rRNA read removal

### Alignment & Processing

-   **STAR** --- genome indexing and read alignment\
-   **SAMtools** --- BAM conversion, sorting, indexing\
-   **StringTie** --- transcript assembly

### Quantification

-   **featureCounts** --- gene-level count matrix\
-   **Salmon** --- transcript quantification\
-   **Kallisto** --- transcript-level pseudoalignment

### Quality Metrics

-   **RSeQC** --- mapping quality and gene body coverage\
-   **Qualimap** --- RNA-seq mapping & coverage metrics

### (Optional) Metagenomics

-   **Kraken2** + **Bracken** --- taxonomic classification

------------------------------------------------------------------------

## 📁 Directory Structure

After a successful run, results are organized under `results/`:

    results/
     ├── fastqc/
     ├── trimmed/
     ├── UMI/
     ├── UMI_dedup/
     ├── aligned/
     ├── bam/
     ├── STAR_index/
     ├── salmon_quant/
     ├── kallisto_index/
     ├── kallisto_quant/
     ├── stringtie/
     ├── featurecounts/
     ├── Qualimap/
     ├── RSeQC/
     ├── multiqc/
     └── ...

------------------------------------------------------------------------

## Install the environment with the tools

### How to install the tools?

    conda env create -f rnaseq.yml

## ⚙️ Pipeline Execution

### Basic Run

    nextflow run main.nf -profile standard

### Use Custom Parameters

    nextflow run main.nf   --reads "/path/to/*_{1,2}.fastq.gz"   --genome "/path/to/genome.fasta"   --gtf "/path/to/annotation.gtf"   --outdir "results"

------------------------------------------------------------------------

## 📦 Software Requirements

The pipeline assumes tools are available in your environment:

  Tool                Purpose
  ------------------- -------------------------------------
  FastQC              Read quality control
  Fastp               Trimming
  STAR                Alignment
  Salmon              Quantification
  Kallisto            Pseudoalignment
  SAMtools            BAM operations
  featureCounts       Gene-level quantification
  StringTie           Transcript assembly
  RSeQC               RNA-seq QC
  Qualimap            Mapping QC
  MultiQC             Summary reports
  SortMeRNA           rRNA removal (optional)
  Kraken2 / Bracken   Taxonomic classification (optional)
  UMI-tools           UMI extraction & deduplication

------------------------------------------------------------------------

## 🧬 Input Requirements

-   Paired-end FASTQ files following this pattern:

```{=html}
<!-- -->
```
    sample1_1.fastq.gz  
    sample1_2.fastq.gz  
    sample2_1.fastq.gz  
    sample2_2.fastq.gz  
    ...

-   Reference genome in FASTA format\
-   GTF annotation file

------------------------------------------------------------------------

## 📤 Output Summary

### Main Outputs

-   Gene counts (`*_counts.txt`)
-   Salmon quantification (`salmon_quant/`)
-   Kallisto quantification (`*_kallisto_output/`)
-   Assembled transcripts (`*_assembled.gtf`)
-   STAR alignment files (SAM/BAM)
-   QC metrics from Qualimap and RSeQC
-   A **MultiQC report** summarizing all results

------------------------------------------------------------------------

## 📚 Citation

If you use this pipeline, please cite the respective tools and Nextflow:

**Di Tommaso P., et al. Nextflow enables reproducible computational
workflows. Nat Biotechnol. 2017.**
