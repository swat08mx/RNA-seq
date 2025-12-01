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
params.reads = "/media/molmed/Analysis/swattik/RNAseq/data/*_{1,2}.fastq.gz"
params.genome = "/media/molmed/Analysis/swattik/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
params.gtf = "/media/molmed/Analysis/swattik/gencode.v49.chr_patch_hapl_scaff.annotation.gtf"
params.outdir = "results"
params.star_index = "/media/molmed/Analysis/swattik/RNAseq/RNA-seq/star_index/"
params.salmon_quant_index = "salmon_quant_index"
//params.kraken_db = "-------------------------------------------------------------------insert-----------------------------------------------------------"
params.qualimap_outdir = "qualimap_outdir"

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
    fastp -i ${sample_id}_1.fastq.gz -o ${sample_id}_trimmed_1.fastq.gz -I ${sample_id}_2.fastq.gz -O ${sample_id}_trimmed_2.fastq.gz -V
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
    path star_index
    tuple val(sample_id), path(reads)
        
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

    script:
    """
    umi_tools extract \\
    --bc-pattern=CCCCCCCCNNNNNNNN \\
    --stdin=${reads[0]} \\
    --stdout=${sample_id}_1_extracted.fastq.gz \\
    --read2-in=${reads[1]} \\
    --read2-out=${sample_id}_2_extracted.fastq.gz \\
    """
}

process DEDUP {
    tag "Deduping ${sample_id}"
    publishDir "${params.outdir}/Dedup", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: dedup_bam
    path("${sample_id}.marked_dup_metrics.txt")

    script:
    """
    gatk MarkDuplicates -I ${sorted_bam} -O ${sample_id}.dedup.bam -M ${sample_id}.marked_dup_metrics.txt
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
    tuple val(sample_id), path("${sorted_bam}.bai"), emit: indexed_file
    
    script:
    """
    samtools index ${sorted_bam}
    """
}

process RSEQC {
    tag "RSEQC ${sample_id}"
    publishDir "${params.outdir}/RSeQC", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_bam_stats.txt")

    script:
    """
    bam_stat.py -i ${bam} > ${sample_id}_bam_stats.txt
    """
}

process QUALIMAP {
    tag "QUALIMAP on ${sample_id}"
    publishDir "${params.outdir}/Qualimap", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path gtf

    output:
    path "${sample_id}_qualimap_outdir"

    script:
    """
    qualimap rnaseq -bam ${bam} -gtf ${gtf} -outdir ${sample_id}_qualimap_outdir -outformat HTML
    """
}

process KRAKEN {
    tag "KRAKEN on ${sample_id}"
    publishDir "${params.outdir}/Kraken", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path kraken_db

    output:
    path "${sample_id}_kraken.out"
    path "${sample_id}.report", emit: kraken_report

    script:
    """
    kraken2 --db $kraken_db --threads 8 --paired --output ${sample_id}_kraken.out \\
    --report ${sample_id}.report ${reads[0]} ${reads[1]}
    """
}

process BRACKEN {
    tag "BRACKEN on ${sample_id}"
    publishDir "${params.outdir}/Bracken", mode: 'copy'

    input:
    path kraken_db
    tuple val(sample_id), path(kraken_report)

    output:
    path "${sample_id}_species_abundance.bracken"

    script:
    """
    bracken -d ${kraken_db} -i ${kraken_report} -o ${sample_id}_species_abundance.bracken -r 100 -l S

    """
}

process KALLISTO_INDEX {
    tag "Kallisto index on ${sample_id}"
    publishDir "${params.outdir}/Kallisto_index", mode: 'copy'

    input:
    path(genome)

    output:
    path "transcripts.idx", emit: transcript_idx

    script:
    """
    kallisto index -i transcripts.idx ${genome}

    """
}

process KALLISTO_QUANT {
    tag "Kallisto quant on ${sample_id}"
    publishDir "${params.outdir}/Kallisto_quant", mode: 'copy'

    input:
    path(transcripts)
    path(reads)

    output:
    path "${sample_id}_kallisto_output/*"

    script:
    """
    kallisto quant -i ${transcripts} -o ${sample_id}_kallisto_output -t 8 ${reads[0]} ${reads[1]}

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
    path gtf
    tuple val(sample_id), path(bam)
        
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
    star_index_ch = Channel.fromPath(params.star_index, checkIfExists: true)
    //kraken_db_ch = Channel.fromPath(params.kraken_db, checkIfExists: true)
    
    // FastQC on raw reads
    FASTQC(reads_ch)
    
    // Trim reads
    FASTP(reads_ch)
    
    // UMI extraction
    UMI_EXTRACT(FASTP.out.trimmed_reads)    
    
    // Salmon index
    SALMON_INDEX(genome_ch)

    // Kallisto index
    //KALLISTO_INDEX(genome_ch)

    // Salmon quantify
    SALMON(UMI_EXTRACT.out.UMI_extracted_reads, SALMON_INDEX.out.salmon_quant_index)

    // Build STAR index
    //BUILD_STAR_INDEX(genome_ch, gtf_ch)
    
    // Align trimmed reads
    STAR_ALIGN(star_index_ch, UMI_EXTRACT.out.UMI_extracted_reads)

    // Convert SAM to BAM
    SAM_TO_BAM(STAR_ALIGN.out.sam)

    // Sort BAM
    SORT_BAM(SAM_TO_BAM.out.bam)

    // Index BAM
    INDEX_BAM(SORT_BAM.out.sorted_bam)

    // Dedup
    DEDUP(SORT_BAM.out.sorted_bam)

    // RSeQC 
    RSEQC(DEDUP.out.dedup_bam)

    // Qualimap
    QUALIMAP(DEDUP.out.dedup_bam, gtf_ch)

    // Kallisto quant
    //KALLISTO_QUANT(KALLISTO_INDEX.out.transcript_idx, reads_ch)
    
    // Kracken
    // KRAKEN(reads_ch, kraken_db_ch)

    // Bracken
    // BRACKEN(kraken_db_ch, KRAKEN.out.kraken_report)

    // Stringtie quant
    STRINGTIE(DEDUP.out.dedup_bam, gtf_ch)

    // Count features
    FEATURECOUNTS(gtf_ch, DEDUP.out.dedup_bam)

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
