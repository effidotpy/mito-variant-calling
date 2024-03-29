process ALIGN_BWA {
    container 'mapping_gatk'

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
    container 'mapping_gatk'

    input:
    tuple val(meta), path(unaligned), path(aligned)
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
    container 'mapping_gatk'

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
        ch_to_merge = ch_revert_bam.join(ch_bam, failOnMismatch:true)
//         ch_to_merge.view()
        ch_merged = MERGE_BAMS(ch_to_merge, chrM_genome)
        ch_sorted = SORT_BAM(ch_merged)

    emit:
        sorted = ch_sorted  // channel: [ val(meta), path(bam) ]
}
