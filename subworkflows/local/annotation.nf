
process HAPLOCHECK {
    container 'mapping_gatk'

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
    container 'mapping_gatk'

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
    container 'mapping_gatk'

    input:
    tuple val(meta), path(vcf)
    val (reference_gnomad)
    val (annoc_file)

    output:
    tuple val(meta), path("annotated_gnomad.vcf")

    script:
    """
    mkdir -p null/temp
    java -jar /app/VarNote-1.2.0.jar Annotation -Q ${vcf} -D:db,tag=GnomadMitochondrial,mode=1 ${reference_gnomad} \\
    -A ${annoc_file} -O annotated_gnomad.vcf -Z false
    """
}


workflow ANNOTATION {

    take:
        ch_vcf          // channel: [ val(meta), path(vcf) ]
        val_gnomad      // string: path to gnomAD VCF database in docker image
        val_annoc       // string: path to VarNote annoc file in docker image

    main:
        ch_gnomad = VARNOTE(ch_vcf, val_gnomad, val_annoc)
        ch_haplocheck = HAPLOCHECK(ch_vcf)
        ch_haplogrep = HAPLOGREP(ch_vcf)

    emit:
        gnomad_annotated = ch_gnomad    // channel: [ val(meta), path(vcf) ]
        haplocheck = ch_haplocheck      // channel: [ val(meta), path(tsv) ]
        haplogrep = ch_haplogrep        // channel: [ val(meta), path(tsv) ]
}