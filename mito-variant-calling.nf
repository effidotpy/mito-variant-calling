params.samplesheet = "$projectDir/data/samplesheet_dev.tsv"
params.outdir = "$projectDir/results"
params.chrM_fasta = "/references/Homo_sapiens_assembly38.chrM.fasta"
params.chrM_shift_fasta = "/references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"


log.info """\
    MITOCHONDRIAL VARIANT CALLING - NF
    ==================================
    samplesheet: ${params.samplesheet}
    """
    .stripIndent()


process NAME_SORT {
    tag "Sorting reads in $bam_file by name"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.sample_id}_namesort.bam")

    script:
    """
    java -jar /app/picard.jar SortSam -I ${bam_file} -O ${meta.sample_id}_namesort.bam -SO queryname --CREATE_INDEX true
    """
}

process FIX_PAIRING {
    tag "Fixing mate pair information in $bam_file"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.sample_id}_fixmate.bam")

    script:
    """
    samtools fixmate ${bam_file} ${meta.sample_id}_fixmate.bam
    """


}

process REMOVE_SINGLETONS {
    tag "Removing singletons from $bam_file"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.sample_id}_singlt.bam")

    script:
    """
    samtools view -h -f 1 -q 1 -O BAM -o ${meta.sample_id}_singlt.bam ${bam_file}
    """
}

process REVERT_BAM {
    //TODO adjust --REMOVE_DUPLICATE_INFORMATION parameter if new samples are included
    tag "Removing alignment information from $bam_file"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.sample_id}_revert.bam")

    script:
    """
    java -jar /app/picard.jar RevertSam -I ${bam_file} -O ${meta.sample_id}_revert.bam \
    --ATTRIBUTE_TO_CLEAR FT --ATTRIBUTE_TO_CLEAR CO --SORT_ORDER queryname --OUTPUT_BY_READGROUP false \
    --RESTORE_ORIGINAL_QUALITIES false --REMOVE_DUPLICATE_INFORMATION false
    """
}


process BAM_TO_FASTQ {
    //TODO adjust --REMOVE_DUPLICATE_INFORMATION parameter if new samples are included
    tag "Extracting reads from $bam_file to (interleaved) FASTQ format"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.sample_id}_mito.fastq")

    script:
    """
    java -jar /app/picard.jar SamToFastq -I ${bam_file} -F ${meta.sample_id}_mito.fastq \
    --INTERLEAVE true --INCLUDE_NON_PF_READS true
    """
}

process ALIGN {
    tag "Aligning $fastq_file to chrM"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(fastq_file)
    val (reference_fasta)

    output:
    tuple val(meta), path("${meta.sample_id}_chrM.bam")

    script:
    """
    ls ${reference_fasta};
    bwa mem -p -Y ${reference_fasta} ${fastq_file} > ${meta.sample_id}_chrM.bam
    """
}

process ALIGN_SHIFT {
    tag "Aligning $fastq_file to chrM"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(fastq_file)
    val (reference_fasta)

    output:
    tuple val(meta), path("${meta.sample_id}_chrM_shift.bam")

    script:
    """
    ls ${reference_fasta};
    bwa mem -p -Y ${reference_fasta} ${fastq_file} > ${meta.sample_id}_chrM_shift.bam
    """
}

process MERGE_BAMS {
    tag "Merging unaligned $unal_bam to with aligned $al_bam"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(unal_bam)
    tuple val(meta), path(al_bam)
    val (reference_fasta)

    output:
    tuple val(meta), path("${meta.sample_id}_merged.bam")

    script:
    """
    java -jar /app/picard.jar MergeBamAlignment -UNMAPPED ${unal_bam} -ALIGNED ${al_bam} -R ${reference_fasta} \\
    -O ${meta.sample_id}_merged.bam --ATTRIBUTES_TO_RETAIN "X0" --ATTRIBUTES_TO_REMOVE NM --ATTRIBUTES_TO_REMOVE MD
    """
}

process MERGE_BAMS_SHIFT {
    tag "Merging unaligned $unal_bam to with aligned $al_bam"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(unal_bam)
    tuple val(meta), path(al_bam)
    val (reference_fasta)

    output:
    tuple val(meta), path("${meta.sample_id}_merged_shift.bam")

    script:
    """
    java -jar /app/picard.jar MergeBamAlignment -UNMAPPED ${unal_bam} -ALIGNED ${al_bam} -R ${reference_fasta} \\
    -O ${meta.sample_id}_merged_shift.bam --ATTRIBUTES_TO_RETAIN "X0" --ATTRIBUTES_TO_REMOVE NM --ATTRIBUTES_TO_REMOVE MD
    """
}

process SORT_BAM {
    tag "Sorting and indexing $bam_file"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.sample_id}_sort.bam")

    script:
    """
    java -jar /app/picard.jar SortSam -I ${bam_file} -O ${meta.sample_id}_sort.bam -SO coordinate --CREATE_INDEX true
    """
}

process SORT_BAM_SHIFT {
    tag "Sorting and indexing $bam_file"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.sample_id}_shift_sort.bam")

    script:
    """
    java -jar /app/picard.jar SortSam -I ${bam_file} -O ${meta.sample_id}_shift_sort.bam -SO coordinate --CREATE_INDEX true
    """
}

workflow {
    Channel
        .fromPath(params.samplesheet)
        | splitCsv( header: true , sep: "\t")
        | map { row ->
            meta = row.subMap('sample_id', 'replicate', 'family' )
            [meta, file(row.bam_file, checkIfExists: true)]
        }
        | set { samples }

    sort_bam_ch = NAME_SORT(samples)
    fix_pair_ch = FIX_PAIRING(sort_bam_ch)
    rmv_singlt_ch = REMOVE_SINGLETONS(fix_pair_ch)
    revert_bam_ch = REVERT_BAM(rmv_singlt_ch)
    mito_fastq_ch = BAM_TO_FASTQ(revert_bam_ch)
    bam_chrM_ch = ALIGN(mito_fastq_ch, params.chrM_fasta)
    bam_chrM_shift_ch = ALIGN_SHIFT(mito_fastq_ch, params.chrM_shift_fasta)
    merged_chrM_ch = MERGE_BAMS(revert_bam_ch, bam_chrM_ch, params.chrM_fasta)
    merged_chrM_shift_ch = MERGE_BAMS_SHIFT(revert_bam_ch, bam_chrM_shift_ch, params.chrM_shift_fasta)
    sorted_ch = SORT_BAM(merged_chrM_ch)
    sorted_shift_ch = SORT_BAM_SHIFT(merged_chrM_shift_ch)
}

