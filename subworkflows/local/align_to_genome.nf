

process ALIGN_BWA {

    input:
    tuple val(meta), path(fastq)
    val (reference_fasta)

    output:
    tuple val(meta), path("aligned.bam")

    script:
    """
    bwa mem -p -Y ${reference_fasta} ${fastq} > aligned.bam
    """
}


process MERGE_BAMS {

    input:
    tuple val(meta), path(unaligned)
    tuple val(meta), path(aligned)
    val (reference_fasta)

    output:
    tuple val(meta), path("merged.bam")

    script:
    """
    java -jar /app/picard.jar MergeBamAlignment -UNMAPPED ${unaligned} -ALIGNED ${aligned} -R ${reference_fasta} \\
    -O merged.bam --ATTRIBUTES_TO_RETAIN "X0" --ATTRIBUTES_TO_REMOVE NM --ATTRIBUTES_TO_REMOVE MD
    """
}

process SORT_BAM {

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("sorted.bam")

    script:
    """
    java -jar /app/picard.jar SortSam -I ${bam} -O sorted.bam -SO coordinate --CREATE_INDEX true
    """
}


workflow ALIGN_TO_GENOME {

    take:
        ch_fastq
        ch_revert_bam
        chrM_genome

    main:
        ch_bam = ALIGN_BWA(ch_fastq, chrM_genome)
        ch_merged = MERGE_BAMS(ch_revert_bam, ch_bam, chrM_genome)
        ch_sorted = SORT_BAM(ch_merged)

    emit:
        sorted = ch_sorted  // channel: [ val(meta), path(bam) ]
}
