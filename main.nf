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
params.gtf = "/home/user1/RNA-seq/rnaseqpipeline/gencode.v49.chr_patch_hapl_scaff.annotation.gtf"
params.outdir = "results"
params.star_index = "star_index"
params.salmon_quant_index = "salmon_quant_index"

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

process SORTMERNA {
    tag "Sorting Rb RNA ${sample_id}"
    publishDir "${params.outdir}/SORTMERNA", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "aligned_out"
    path "other_out"

    script:
    """
    sortmerna --ref ${params.genome} --reads ${reads[0]} \\
    --reads ${reads[1]} --aligned ${sample_id}_rRNA_reads --other ${sample_id}_non_rRNA_reads --fastx --paired_out

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
         --sjdbOverhang 99 \\
         --genomeSAsparseD 2 \\
         --limitGenomeGenerateRAM 10781458698
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
        --genomeSAsparseD 2 \\
        --outFileNamePrefix ${sample_id}_
    """
}

process SALMON_INDEX {
    tag "Salmon index"
    publishDir "${params.outdir}/${params.salmon_quant_index}", mode: 'copy', pattern: "*.log"

    input:
    path genome

    output:
    path "${params.salmon_quant_index}", emit: salmon_quant_index

    script:
    """
    salmon index -t ${genome} -i ${params.salmon_quant_index}
    """
}


process SALMON {
    tag "Salmon quantification ${sample_id}"
    publishDir "${params.outdir}/salmon_quant", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(reads)
    path salmon_quant_index

    output:
    path "salmon_quant", emit: salmon_quant

    script:
    """
    salmon quant -i ${params.salmon_quant_index} -l A -1 ${reads[0]} -2 ${reads[1]} \\
    --validateMappings -o salmon_quant
    """
}

process UMI_EXTRACT {
    tag "UMI extraction on ${sample_id}"
    publishDir "${params.outdir}/UMI", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_*_extracted.fastq.gz"), emit: UMI_extracted_reads
    path "${sample_id}_whitelist.txt", emit: whitelist

    script:
    """
    umi_tools extract \\
    --bc-pattern=CCCCCCCCNNNNNNNN \\
    --stdin=${reads[0]} \\
    --stdout=${sample_id}_1_extracted.fastq.gz \\
    --read2-in=${reads[1]} \\
    --read2-out=${sample_id}_2_extracted.fastq.gz \\
    --filter-cell-barcode \\
    --whitelist=${sample_id}_whitelist.txt
    """
}

process UMI_DEDUP {
    tag "UMI deduping ${sample_id}"
    publishDir "${params.outdir}/UMI_dedup", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: dedup_bam

    script:
    """
    umi_tools dedup \\
    -I ${sorted_bam} \\
    -S ${sample_id}.dedup.bam \\
    --paired
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

process STRINGTIE {
    tag "Stringtie on ${sample_id}"
    publishDir "${params.outdir}/stringtie", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}_assembled.gtf"), emit: assembled_gtf

    script:
    """
    stringtie -p 8 -G ${gtf} -o ${sample_id}_assembled.gtf ${sorted_bam}

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
    
    // UMI extraction
    UMI_EXTRACT(FASTP.out.trimmed_reads)    
    
    // Salmon index
    SALMON_INDEX(genome_ch)

    // Salmon quantify
    SALMON(FASTP.out.trimmed_reads, SALMON_INDEX.out.salmon_quant_index)

    // Build STAR index
    BUILD_STAR_INDEX(genome_ch, gtf_ch)
    
    // Align trimmed reads
    STAR_ALIGN(FASTP.out.trimmed_reads, BUILD_STAR_INDEX.out.index)

    // Convert SAM to BAM
    SAM_TO_BAM(STAR_ALIGN.out.sam)

    // Sort BAM
    SORT_BAM(SAM_TO_BAM.out.bam)

    // Dedup and UMI removal
    UMI_DEDUP(SORT_BAM.out.sorted_bam)

    // Index BAM
    INDEX_BAM(UMI_DEDUP.out.dedup_bam)

    // Stringtie quant
    STRINGTIE(INDEX_BAM.out.indexed_bam, gtf_ch)

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
