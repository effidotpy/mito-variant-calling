process VAR_CALLING {
    container 'mapping_gatk'

    input:
    tuple val(meta), path(bam)
    val (reference_fasta)

    output:
    tuple val(meta), path("mutect2_ref.vcf")        , emit: vcf
    tuple val(meta), path("mutect2_ref.vcf.stats")  , emit: stats

    script:
    """
    gatk Mutect2 -I ${bam} -O mutect2_ref.vcf --reference ${reference_fasta} -L chrM:576-16024 \\
    --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter \\
    --annotation StrandBiasBySample --max-reads-per-alignment-start 0 --max-mnp-distance 0 --mitochondria-mode
    """
}

process VAR_CALLING_SHIFT {
    container 'mapping_gatk'

    input:
    tuple val(meta), path(bam)
    val (reference_fasta)

    output:
    tuple val(meta), path("mutect2_shift.vcf")          , emit: vcf
    tuple val(meta), path("mutect2_shift.vcf.stats")    , emit: stats

    script:
    """
    gatk Mutect2 -I ${bam} -O mutect2_shift.vcf --reference ${reference_fasta} -L chrM:8025-9144 \\
    --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter \\
    --annotation StrandBiasBySample --max-reads-per-alignment-start 0 --max-mnp-distance 0 --mitochondria-mode
    """
}

process LIFTOVER {
    container 'mapping_gatk'

    input:
    tuple val(meta), path(vcf)
    val (reference_fasta)
    val (shiftback_chain)

    output:
    tuple val(meta), path("shiftback.vcf")

    script:
    """
    java -jar /app/picard.jar LiftoverVcf -I ${vcf} -O shiftback.vcf -C ${shiftback_chain} -R ${reference_fasta} \\
    --REJECT rejected.vcf
    """
}

process MERGE_VCFS {
    container 'mapping_gatk'

    input:
    tuple val(meta), path(vcf_ref), path(vcf_shift)

    output:
    tuple val(meta), path("merged.vcf")

    script:
    """
    java -jar /app/picard.jar MergeVcfs -I ${vcf_ref} -I ${vcf_shift} -O merged.vcf
    """
}

process MERGE_STATS {
    container 'mapping_gatk'

    input:
    tuple val(meta), path (stats_ref), path (stats_shift)

    output:
    tuple val(meta), path("merged.vcf.stats")

    script:
    """
    gatk MergeMutectStats --stats ${stats_ref} --stats ${stats_shift} -O merged.vcf.stats
    """
}

process FILTER_MUTECT {
    container 'mapping_gatk'

    input:
    tuple val(meta), path(vcf), path(stats)
    val (reference_fasta)

    output:
    tuple val(meta), path("filtered.vcf")

    script:
    """
    gatk FilterMutectCalls -V ${vcf} --stats ${stats} -O filtered.vcf -R ${reference_fasta} \\
    --mitochondria-mode
    """
}

process BLACKLIST_SITES {
    container 'mapping_gatk'

    input:
    tuple val(meta), path(vcf)
    val (blacklist)

    output:
    tuple val(meta), path("blacklisted.vcf")

    script:
    """
    gatk VariantFiltration -V ${vcf} -O blacklisted.vcf --mask ${blacklist} --mask-name blacklisted
    """
}

process NORMALIZE_CALLS {
    container 'mapping_gatk'

    input:
    tuple val(meta), path(vcf)
    val (reference_fasta)

    output:
    tuple val(meta), path("normalized.vcf")

    script:
    """
    gatk LeftAlignAndTrimVariants -V ${vcf} -R ${reference_fasta} --split-multi-allelics --dont-trim-alleles \\
    --keep-original-ac -O normalized.vcf
    """
}

workflow VARIANT_CALLING {

    take:
        ch_bam_ref          // channel: [ val(meta), path(bam) ]
        ch_bam_shift        // channel: [ val(meta), path(bam) ]
        val_genome_ref      // string:  path to reference chrM genome
        val_genome_shift    // string:  path to (shift) reference chrM genome
        val_shiftback       // string:  path to shiftback chain file
        val_blacklist       // string:  path to variant's black list


    main:
        ch_var_ref = VAR_CALLING(ch_bam_ref, val_genome_ref)
        ch_var_shift = VAR_CALLING_SHIFT(ch_bam_shift, val_genome_shift)
        ch_shiftback = LIFTOVER(ch_var_shift.vcf, val_genome_ref, val_shiftback)

        ch_vcfs_to_merge = ch_var_ref.vcf.join(ch_shiftback, failOnMismatch:true)
        ch_vcf_merged = MERGE_VCFS(ch_vcfs_to_merge)

        ch_stats_to_merge = ch_var_ref.stats.join(ch_var_shift.stats, failOnMismatch:true)
        ch_stats_merged = MERGE_STATS(ch_stats_to_merge)

        ch_to_filter = ch_vcf_merged.join(ch_stats_merged, failOnMismatch:true)
        ch_filter_mutect = FILTER_MUTECT(ch_to_filter, val_genome_ref)

        ch_blacklisted = BLACKLIST_SITES(ch_filter_mutect, val_blacklist)
        ch_normalized = NORMALIZE_CALLS(ch_blacklisted, val_genome_ref)

    emit:
        normalized = ch_normalized   // channel: [ val(meta), path(vcf) ]




}