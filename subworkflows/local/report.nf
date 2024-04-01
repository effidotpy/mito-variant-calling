
process VCF2TSV {
    conda "reporting.yml"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    path("*.tsv")

    script:
    """
    vcf2tsv.py --vcf_in ${vcf} --tsv_out ${meta.sample_id}_variants.tsv
    """
}

process GENERATE_HTML {
    conda "reporting.yml"
    publishDir params.outdir, mode: 'copy'

    input:
    path tsv_variants
    path tsv_haplocheck
    path tsv_haplogrep
    path tsv_depth

    output:
    path("report.html")

    script:
    """
    generate_report.py --var ${tsv_variants} --hc ${tsv_haplocheck} --hg ${tsv_haplogrep} --dp ${tsv_depth} --out report.html
    """
}

workflow REPORT {

    take:
    ch_vcf
    ch_haplocheck
    ch_haplogrep
    ch_depth

    main:
    ch_tsv = VCF2TSV(ch_vcf)
    GENERATE_HTML ( ch_tsv.collect(), ch_haplocheck.collect(), ch_haplogrep.collect(), ch_depth.collect() )

    emit:
    tsv = ch_tsv

}