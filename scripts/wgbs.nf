#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ───────────── PARAMETERS ─────────────
params.fastq         = "${workflow.projectDir}/data/*_{1,2}.fastq"  // Matches paired-end files
params.results       = "${workflow.projectDir}/results"

params.adapter       = 'AGATCGGAAGAGC'
params.quality       = 20
params.min_length    = 50
params.threads       = 4
params.genome        = "${workflow.projectDir}/genome"  // Path to Bismark genome index

// ───────────── PROCESS: GENOME PREPARATION ─────────────
process GENOME_PREP {
    tag "Genome preparation"

    input:
    path genome_dir

    output:
    path genome_dir, emit: prepared_genome

    script:
    """
    # Check if Bismark genome is already prepared
    if [ -d "${genome_dir}/Bisulfite_Genome" ]; then
        echo "Bismark genome already prepared in ${genome_dir}"
    else
        echo "Preparing Bismark genome index..."
        bismark_genome_preparation --bowtie2 ${genome_dir}
        echo "Genome preparation completed"
    fi
    """
}

// ───────────── PROCESS: FASTQC on raw reads ─────────────
process QC {
    tag "$sample_id"
    publishDir "${params.results}/${sample_id}/01_QC_reports", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc -t ${params.threads} ${read1} ${read2}
    """
}

// ───────────── PROCESS: CUTADAPT ─────────────
process TRIM {
    tag "$sample_id"
    publishDir "${params.results}/${sample_id}/02_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_cutadapt.log", emit: log

    script:
    """
    cutadapt \\
        -j ${params.threads} \\
        -a ${params.adapter} \\
        -A ${params.adapter} \\
        -q ${params.quality} \\
        --minimum-length ${params.min_length} \\
        -o ${sample_id}_R1.trimmed.fastq.gz \\
        -p ${sample_id}_R2.trimmed.fastq.gz \\
        ${read1} ${read2} > ${sample_id}_cutadapt.log 2>&1
    """
}

// ───────────── PROCESS: FASTQC on trimmed reads ─────────────
process TRIM_QC {
    tag "$sample_id"
    publishDir "${params.results}/${sample_id}/03_trimmed_QC_reports", mode: 'copy'

    input:
    tuple val(sample_id), path(read1_trimmed), path(read2_trimmed)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc -t ${params.threads} ${read1_trimmed} ${read2_trimmed}
    """
}

// ───────────── PROCESS: BISMARK ALIGNMENT ─────────────
process ALIGN {
    tag "$sample_id"
    publishDir "${params.results}/${sample_id}/04_aligned", mode: 'copy'
    
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(sample_id), path(read1_trimmed), path(read2_trimmed), path(prepared_genome)

    output:
    tuple val(sample_id), path("${sample_id}_aligned.bam"), emit: bam
    path "*.txt", emit: report

    script:
    """
    echo "Starting alignment for ${sample_id}"
    echo "Input files: ${read1_trimmed}, ${read2_trimmed}"
    echo "Genome: ${prepared_genome}"
    
    bismark \\
        --multicore \$((${params.threads}/2)) \\
        --genome ${prepared_genome} \\
        -1 ${read1_trimmed} \\
        -2 ${read2_trimmed} \\
        -o .

    echo "Bismark alignment completed for ${sample_id}"
    echo "Looking for BAM files..."
    ls -la *.bam
    
    # Find and rename the BAM file
    bam_file=\$(find . -name "*.bam" -type f | head -1)
    if [ -n "\$bam_file" ]; then
        echo "Found BAM file: \$bam_file"
        mv "\$bam_file" "${sample_id}_aligned.bam"
        echo "Renamed to: ${sample_id}_aligned.bam"
    else
        echo "ERROR: No BAM file found after alignment!"
        echo "Files in directory:"
        ls -la
        exit 1
    fi
    """
}

// ───────────── PROCESS: DEDUPLICATION ─────────────
process DEDUP {
    tag "$sample_id"
    publishDir "${params.results}/${sample_id}/05_deduplicated", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_deduplicated.bam"), emit: dedup_bam
    path "*.txt", emit: report

    script:
    """
    deduplicate_bismark -p --bam ${bam}
    
    # Find and rename the deduplicated BAM file
    dedup_file=\$(find . -name "*deduplicated.bam" -type f | head -1)
    if [ -n "\$dedup_file" ]; then
        mv "\$dedup_file" "${sample_id}_deduplicated.bam"
    else
        echo "Error: No deduplicated BAM file found"
        exit 1
    fi
    """
}

// ───────────── PROCESS: METHYLATION EXTRACTION ─────────────
process METH_EXTRACT {
    tag "$sample_id"
    publishDir "${params.results}/${sample_id}/06_methylation", mode: 'copy'

    input:
    tuple val(sample_id), path(dedup_bam), path(prepared_genome)

    output:
    path "*.txt.gz"
    path "*.bedGraph.gz", optional: true
    path "*.cov.gz", optional: true
    path "*_splitting_report.txt", optional: true
    path "*_mbias.txt", optional: true
    path "*.cytosine_context_summary.txt", optional: true
    path "*CpG_report.txt.gz", optional: true

    script:
    """
    bismark_methylation_extractor \\
        --paired-end \\
        --gzip \\
        --bedGraph \\
        --cytosine_report \\
        --comprehensive \\
        --merge_non_CpG \\
        --genome_folder ${prepared_genome} \\
        ${dedup_bam} \\
        -o .
    """
}

// ───────────── WORKFLOW ─────────────
workflow {

    log.info "🔍 Looking for FASTQ files matching: ${params.fastq}"

    // Prepare genome index if needed
    genome_channel = Channel.fromPath(params.genome, checkIfExists: true, type: 'dir')
    prepared_genome = GENOME_PREP(genome_channel)

    // Load paired-end FASTQ files
    fastq_pairs = Channel.fromFilePairs(params.fastq, checkIfExists: true)
        .map { sample_id, reads -> 
            log.info "📂 Loaded: ${sample_id}, ${reads[0].getName()}, ${reads[1].getName()}"
            tuple(sample_id, reads[0], reads[1])
        }

    // Run FastQC on raw reads
    fastq_pairs | QC

    // Run Cutadapt
    TRIM(fastq_pairs)

    // Run FastQC on trimmed reads
    TRIM.out.trimmed_reads | TRIM_QC

    // Combine trimmed reads with prepared genome for each sample
    trimmed_with_genome = TRIM.out.trimmed_reads
        .combine(prepared_genome.prepared_genome)

    // Run Bismark alignment
    ALIGN(trimmed_with_genome)

    // Run deduplication
    DEDUP(ALIGN.out.bam)

    // Combine dedup output with prepared genome for methylation extraction
    dedup_with_genome = DEDUP.out.dedup_bam
        .combine(prepared_genome.prepared_genome)

    // Run methylation extraction
    METH_EXTRACT(dedup_with_genome)
}