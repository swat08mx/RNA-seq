#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    RNA-seq Pipeline with featureCounts
========================================================================================
    A containerless RNA-seq pipeline including:
    - FastQC for quality control
    - Fastp for adapter trimming
    - STAR for alignment
    - SAMtools for BAM processing
    - featureCounts for read quantification
========================================================================================
*/

// Parameters
params.reads = "/home/user1/RNA-seq/data/*_{1,2}.fastq.gz"
params.genome = "/home/user1/test/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
params.gtf = "/home/user1/RNA-seq/rnaseq-pipeline/gencode.v49.chr_patch_hapl_scaff.annotation.gtf"
params.outdir = "results"
params.star_index = "star_index"

// Print parameter summary
log.info """\
    R N A - S E Q   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    genome       : ${params.genome}
    gtf          : ${params.gtf}
    star_index : ${params.star_index}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process FASTQC {
    tag "FastQC on ${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "*.{html,zip}"
    
    script:
    """
    fastqc -q ${reads}
    """
}

process FASTP {
    tag "Trimming ${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*.fastq.gz"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_*.fastq.gz"), emit: trimmed_reads
    
    script:
    """
    /home/user1/test/tools/fastp -i ${sample_id}_1.fastq.gz -o ${sample_id}_trimmed_1.fastq.gz -I ${sample_id}_2.fastq.gz -O ${sample_id}_trimmed_2.fastq.gz -V
    """
}


process BUILD_STAR_INDEX {
    tag "Building STAR index"
    publishDir "${params.outdir}/${params.star_index}", mode: 'copy'
    
    input:
    path genome
    path gtf

    output:
    path "${params.star_index}", emit: index
    
      
    script:
    """
    mkdir -p ${params.star_index}
    STAR --runThreadN ${task.cpus} \\
         --runMode genomeGenerate \\
         --genomeDir ${params.star_index} \\
         --genomeFastaFiles ${genome} \\
         --sjdbGTFfile ${gtf} \\
         --sjdbOverhang 99
    """
}

process STAR_ALIGN {
    tag "Aligning ${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy', pattern: "*.log"
    
    input:
    tuple val(sample_id), path(reads)
    path star_index
    
    output:
    tuple val(sample_id), path("${sample_id}_Aligned.out.sam"), emit: sam
    path "${sample_id}_Log.final.out", emit: log
    
    script:
    """
    STAR --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample_id}_
    """
}

process SAM_TO_BAM {
    tag "Converting ${sample_id} SAM to BAM"
    publishDir "${params.outdir}/bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    
    script:
    """
    samtools view -@ ${task.cpus} -bS ${sam} > ${sample_id}.bam
    """
}

process SORT_BAM {
    tag "Sorting ${sample_id}"
    publishDir "${params.outdir}/bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: sorted_bam
    
    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}_sorted.bam ${bam}
    """
}

process INDEX_BAM {
    tag "Indexing ${sample_id}"
    publishDir "${params.outdir}/bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sorted_bam)
    
    output:
    tuple val(sample_id), path(sorted_bam), path("${sorted_bam}.bai"), emit: indexed_bam
    
    script:
    """
    samtools index ${sorted_bam}
    """
}

process FEATURECOUNTS {
    tag "Counting features for ${sample_id}"
    publishDir "${params.outdir}/featurecounts", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path gtf
    
    output:
    path "${sample_id}_counts.txt", emit: counts
    path "${sample_id}_counts.txt.summary", emit: summary
    
    script:
    """
    featureCounts \\
        -p \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o ${sample_id}_counts.txt \\
        ${bam}
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path '*'
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc .
    """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    // Create input channels
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    genome_ch = params.genome ? Channel.fromPath(params.genome, checkIfExists: true) : Channel.empty()
    gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)

    // FastQC on raw reads
    FASTQC(reads_ch)
    
    // Trim reads
    FASTP(reads_ch)
    
    // Build STAR index
    BUILD_STAR_INDEX(genome_ch, gtf_ch)
    
    // Align trimmed reads
    STAR_ALIGN(FASTP.out.trimmed_reads, BUILD_STAR_INDEX.out.index)
    
    // Convert SAM to BAM
    SAM_TO_BAM(STAR_ALIGN.out.sam)

    // Sort BAM
    SORT_BAM(SAM_TO_BAM.out.bam)

    // Index BAM
    INDEX_BAM(SORT_BAM.out.sorted_bam)

    // Count features
    FEATURECOUNTS(INDEX_BAM.out.indexed_bam, gtf_ch)

    // Collect all QC outputs for MultiQC
    multiqc_input = FASTQC.out
        .mix(STAR_ALIGN.out.log)
        .mix(FEATURECOUNTS.out.summary)
        .collect()

    MULTIQC(multiqc_input)
}


workflow.onComplete {
    log.info """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
}
