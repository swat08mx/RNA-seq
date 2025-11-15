# RNA-seq Pipeline with featureCounts

A comprehensive Nextflow DSL2 RNA-seq analysis pipeline without container dependencies. This pipeline performs quality control, read trimming, alignment, and quantification using featureCounts.

## Pipeline Overview

The pipeline includes the following steps:

1. **FastQC** - Quality control of raw reads
2. **Trimmomatic** - Adapter trimming and quality filtering
3. **HISAT2** - Read alignment to reference genome
4. **SAMtools** - SAM to BAM conversion, sorting, and indexing
5. **featureCounts** - Read quantification at gene level
6. **MultiQC** - Aggregate quality control report

## Prerequisites

The following tools must be installed and available in your PATH:

- Nextflow (>= 23.04.0)
- FastQC (>= 0.11.9)
- Trimmomatic (>= 0.39)
- HISAT2 (>= 2.2.0)
- SAMtools (>= 1.10)
- Subread (for featureCounts) (>= 2.0.0)
- MultiQC (>= 1.9)

### Installation Example (Ubuntu/Debian)

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y fastqc trimmomatic hisat2 samtools

# Install Subread (for featureCounts)
wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-Linux-x86_64.tar.gz
tar -xzf subread-2.0.6-Linux-x86_64.tar.gz
export PATH=$PATH:$(pwd)/subread-2.0.6-Linux-x86_64/bin

# Install MultiQC via pip
pip install multiqc
```

## Required Inputs

1. **Paired-end FASTQ files** - Raw sequencing reads
2. **Reference genome** - FASTA format
3. **Gene annotation** - GTF format
4. **Adapter sequences** - For Trimmomatic (usually included with installation)

## Usage

### Quick Start

```bash
# Run with command-line parameters
nextflow run main.nf \\
    --reads "data/*_R{1,2}.fastq.gz" \\
    --genome reference/genome.fa \\
    --gtf annotation/genes.gtf \\
    --outdir results
```

### Using a Parameters File

```bash
# Copy and edit the parameters file
cp params.yaml my_params.yaml
# Edit my_params.yaml with your file paths

# Run with parameters file
nextflow run main.nf -params-file my_params.yaml
```

### Using Pre-built HISAT2 Index

If you have a pre-built HISAT2 index, you can skip the indexing step:

```bash
nextflow run main.nf \\
    --reads "data/*_R{1,2}.fastq.gz" \\
    --gtf annotation/genes.gtf \\
    --hisat2_index reference/hisat2_index/genome \\
    --outdir results
```

## Input File Naming

The pipeline expects paired-end reads with the following naming pattern:
- `sample1_R1.fastq.gz` and `sample1_R2.fastq.gz`
- `sample2_R1.fastq.gz` and `sample2_R2.fastq.gz`

The `{1,2}` in the glob pattern matches R1 and R2 files automatically.

## Output Structure

```
results/
├── fastqc/                 # FastQC reports
├── trimmed/                # Trimmed FASTQ files
├── hisat2_index/           # HISAT2 index (if built)
├── aligned/                # Alignment logs
├── bam/                    # BAM files (sorted and indexed)
├── featurecounts/          # Gene count matrices
│   ├── *_counts.txt        # Count tables
│   └── *_counts.txt.summary # Counting statistics
├── multiqc/                # Aggregated QC report
└── pipeline_info/          # Execution reports
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

## Key Output Files

- **featureCounts output**: `results/featurecounts/*_counts.txt`
  - Tab-delimited file with gene IDs and read counts
  - Column 1: Gene ID
  - Column 2-6: Gene annotation info
  - Column 7+: Read counts

- **MultiQC report**: `results/multiqc/multiqc_report.html`
  - Comprehensive quality control summary

## Resource Configuration

Default resources are defined in `nextflow.config`. Adjust based on your system:

```groovy
process {
    withName: 'HISAT2_ALIGN' {
        cpus   = 8
        memory = '12.GB'
        time   = '6.h'
    }
}
```

## Execution Profiles

### Local execution (default)
```bash
nextflow run main.nf -params-file params.yaml
```

### SLURM cluster
```bash
nextflow run main.nf -params-file params.yaml -profile slurm
```

### AWS Batch
```bash
nextflow run main.nf -params-file params.yaml -profile aws
```

## Troubleshooting

### Issue: "Command not found"
**Solution**: Ensure all required tools are installed and in your PATH

### Issue: "No such file or directory"
**Solution**: Use absolute paths for input files or verify relative paths

### Issue: HISAT2 index building takes too long
**Solution**: Pre-build the index and provide it via `--hisat2_index`

### Issue: Out of memory errors
**Solution**: Adjust memory allocations in `nextflow.config` based on your data size

## Pipeline Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--reads` | Path to paired-end FASTQ files | `data/*_R{1,2}.fastq.gz` |
| `--genome` | Path to reference genome (FASTA) | `null` |
| `--gtf` | Path to gene annotation (GTF) | `null` |
| `--hisat2_index` | Path to pre-built HISAT2 index | `null` |
| `--outdir` | Output directory | `results` |
| `--adapter` | Adapter file for Trimmomatic | `TruSeq3-PE.fa` |

## featureCounts Options

The pipeline uses these featureCounts parameters:
- `-p`: Count fragments (paired-end mode)
- `-T`: Number of threads
- `-a`: Gene annotation file (GTF)

To customize featureCounts behavior, edit the `FEATURECOUNTS` process in `main.nf`.

Common modifications:
```bash
# Count at exon level instead of gene
featureCounts -t exon -g gene_id ...

# Count multimapping reads
featureCounts -M ...

# Strand-specific counting
featureCounts -s 1 ...  # forward strand
featureCounts -s 2 ...  # reverse strand
```

## Citation

If you use this pipeline, please cite the tools:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology.
- **FastQC**: Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data.
- **Trimmomatic**: Bolger, A.M., et al. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data.
- **HISAT2**: Kim, D., et al. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype.
- **SAMtools**: Li, H., et al. (2009). The Sequence Alignment/Map format and SAMtools.
- **featureCounts**: Liao, Y., et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.
- **MultiQC**: Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report.

## License

This pipeline is available under the MIT License.
