
process ADD_READGROUP {

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("readgroup.bam")

    script:
    """
    samtools addreplacerg -r '@RG\tID:${meta.sample_id}_${meta.source}\tSM:${meta.sample_id}_${meta.source}' \\
    ${bam_file} -o readgroup.bam
    """
}

process NAME_SORT {

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("namesort.bam")

    script:
    """
    java -jar /app/picard.jar SortSam -I ${bam_file} -O namesort.bam -SO queryname --CREATE_INDEX true
    """
}

process FIX_PAIRING {

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("fixmate.bam")

    script:
    """
    samtools fixmate ${bam_file} fixmate.bam
    """


}

process REMOVE_SINGLETONS {

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("singleton.bam")

    script:
    """
    samtools view -h -f 1 -q 1 -O BAM -o singleton.bam ${bam_file}
    """
}

process REVERT_BAM {
    //TODO adjust --REMOVE_DUPLICATE_INFORMATION parameter if new samples are included

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("revert.bam")

    script:
    """
    java -jar /app/picard.jar RevertSam -I ${bam_file} -O revert.bam \
    --ATTRIBUTE_TO_CLEAR FT --ATTRIBUTE_TO_CLEAR CO --SORT_ORDER queryname --OUTPUT_BY_READGROUP false \
    --RESTORE_ORIGINAL_QUALITIES false --REMOVE_DUPLICATE_INFORMATION false
    """
}

process BAM_TO_FASTQ {

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("mito.fastq")

    script:
    """
    java -jar /app/picard.jar SamToFastq -I ${bam_file} -F mito.fastq \
    --INTERLEAVE true --INCLUDE_NON_PF_READS true
    """
}

workflow PROCESS_INPUT_BAM {

    take:
        ch_samplesheet  // channel: [ val(meta), path(bam) ]

    main:
        ch_readgroup = ADD_READGROUP(ch_samplesheet)
        ch_sort_bam = NAME_SORT(ch_readgroup)
        ch_fix_pair = FIX_PAIRING(ch_sort_bam)
        ch_rmv_singlt = REMOVE_SINGLETONS(ch_fix_pair)
        ch_revert_bam = REVERT_BAM(ch_rmv_singlt)
        ch_mito_fastq = BAM_TO_FASTQ(ch_revert_bam)

    emit:
        mito_fastq = ch_mito_fastq  // channel: [ val(meta), path(fastq) ]
        revert_bam = ch_revert_bam  // channel: [ val(meta), path(bam) ]
}
