params.samplesheet = "$projectDir/data/samplesheet_dev.tsv"
params.outdir = "$projectDir/results"
params.chrM_fasta = "/references/Homo_sapiens_assembly38.chrM.fasta"
params.chrM_shift_fasta = "/references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
params.shiftback_chain = "/references/ShiftBack.chain"
params.blacklist_sites = "/references/blacklist_sites.hg38.chrM.bed"
params.reference_gnomad = "/references/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
params.annoc_file = "/references/varnote.annoc"


log.info """\
    MITOCHONDRIAL VARIANT CALLING - NF
    ==================================
    samplesheet: ${params.samplesheet}
    """
    .stripIndent()


process ADD_READGROUP {
    tag "Adding sample information to ${bam_file} readgroup"
    //publishDir params.outdir, mode: 'copy'

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
    tag "Sorting reads in $bam_file by name"
    //publishDir params.outdir, mode: 'copy'

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

process VAR_CALLING {
    tag "Calling variants in $bam_file"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)
    val (reference_fasta)

    output:
    tuple val(meta), path("mutect2_raw.vcf")
    path("mutect2_raw.vcf.stats")

    script:
    """
    gatk Mutect2 -I ${bam_file} -O mutect2_raw.vcf --reference ${reference_fasta} -L chrM:576-16024 \\
    --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter \\
    --annotation StrandBiasBySample --max-reads-per-alignment-start 0 --max-mnp-distance 0 --mitochondria-mode
    """
}

process VAR_CALLING_SHIFT {
    tag "Calling variants in $bam_file"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam_file)
    val (reference_fasta)

    output:
    tuple val(meta), path("mutect2_shift_raw.vcf")
    path("mutect2_shift_raw.vcf.stats")

    script:
    """
    gatk Mutect2 -I ${bam_file} -O mutect2_shift_raw.vcf --reference ${reference_fasta} -L chrM:8025-9144 \\
    --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateUnmappedAndUnmappedReadFilter \\
    --annotation StrandBiasBySample --max-reads-per-alignment-start 0 --max-mnp-distance 0 --mitochondria-mode
    """
}

process LIFTOVER {
    tag "Adjusting positions in $vcf_file"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)
    path (vcf_stats)
    val (reference_fasta)
    val (shiftback_chain)

    output:
    tuple val(meta), path("mutect2_shiftback.vcf")

    script:
    """
    java -jar /app/picard.jar LiftoverVcf -I ${vcf_file} -O mutect2_shiftback.vcf -C ${shiftback_chain} -R ${reference_fasta} \\
    --REJECT rejected.vcf
    """
}

process MERGE_VCFS {
    tag "Merging VCFs ${vcf_file} and ${vcf_file_shift}"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)
    path (vcf_stats)
    tuple val(meta_shift), path(vcf_file_shift)

    output:
    tuple val(meta), path("merged.vcf")

    script:
    """
    java -jar /app/picard.jar MergeVcfs -I ${vcf_file} -I ${vcf_file_shift} -O merged.vcf
    """
}

process MERGE_STATS {
    tag "Merging stats ${vcf_stats} and ${vcf_stats_shift}"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)
    path (vcf_stats)
    tuple val(meta_shift), path(vcf_file_shift)
    path (vcf_stats_shift)

    output:
    tuple val(meta), path("merged.vcf.stats")

    script:
    """
    gatk MergeMutectStats --stats ${vcf_stats} --stats ${vcf_stats_shift} -O merged.vcf.stats
    """
}

process FILTER_MUTECT {
    tag "Filtering calls in ${vcf_file}"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)
    tuple val(meta), path(stats_file)
    val (reference_fasta)

    output:
    tuple val(meta), path("filtered.vcf")

    script:
    """
    gatk FilterMutectCalls -V ${vcf_file} --stats ${stats_file} -O filtered.vcf -R ${reference_fasta} \\
    --mitochondria-mode
    """
}

process BLACKLIST_SITES {
    tag "Annotating blacklisted sites in ${vcf_file}"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)
    val (blacklist)

    output:
    tuple val(meta), path("blacklisted.vcf")

    script:
    """
    gatk VariantFiltration -V ${vcf_file} -O blacklisted.vcf --mask ${blacklist} --mask-name blacklisted
    """
}

process NORMALIZE_CALLS {
    tag "Normalizing variants in ${vcf_file}"
    //publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)
    val (reference_fasta)

    output:
    tuple val(meta), path("normalized.vcf")

    script:
    """
    gatk LeftAlignAndTrimVariants -V ${vcf_file} -R ${reference_fasta} --split-multi-allelics --dont-trim-alleles \\
    --keep-original-ac -O normalized.vcf
    """
}

process HAPLOCHECK {
    tag "Checking for contamination in ${vcf_file}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)

    output:
    tuple val(meta), path("haplocheck.tsv")

    script:
    """
    /app/haplocheck --out haplocheck.tsv --raw ${vcf_file}
    """
}

process HAPLOGREP {
    tag "Inferring haplogoup in ${vcf_file}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)

    output:
    tuple val(meta), path("haplogrep.tsv")

    script:
    """
    /app/haplogrep classify --in ${vcf_file} --out haplogrep.tsv --format vcf --extend-report
    """
}




process VARNOTE {
    tag "Annotating gnomAD in ${vcf_file}"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)
    val (reference_gnomad)
    val (annoc_file)

    output:
    tuple val(meta), path("annotated_gnomad.vcf")

    script:
    """
    mkdir -p null/temp
    java -jar /app/VarNote-1.2.0.jar Annotation -Q ${vcf_file} -D:db,tag=GnomadMitochondrial,mode=1 ${reference_gnomad} \\
    -A ${annoc_file} -O annotated_gnomad.vcf -Z false
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

    readgroup_ch = ADD_READGROUP(samples)
    sort_bam_ch = NAME_SORT(readgroup_ch)
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
    var_ch = VAR_CALLING(sorted_ch, params.chrM_fasta)
    var_shift_ch = VAR_CALLING_SHIFT(sorted_shift_ch, params.chrM_shift_fasta)
    var_shiftback_ch = LIFTOVER(var_shift_ch, params.chrM_fasta, params.shiftback_chain)
    merged_vcf_ch = MERGE_VCFS(var_ch, var_shiftback_ch)
    merged_stats_ch = MERGE_STATS(var_ch, var_shift_ch)
    filter_mutect_ch = FILTER_MUTECT(merged_vcf_ch, merged_stats_ch, params.chrM_fasta)
    blacklisted_ch = BLACKLIST_SITES(filter_mutect_ch, params.blacklist_sites)
    normalized_ch = NORMALIZE_CALLS(blacklisted_ch, params.chrM_fasta)
    annot_gnomad_ch = VARNOTE(normalized_ch, params.reference_gnomad, params.annoc_file)
    haplocheck_ch = HAPLOCHECK(normalized_ch)
    haplogrep_ch = HAPLOGREP(normalized_ch)

}

