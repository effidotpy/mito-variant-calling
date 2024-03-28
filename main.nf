#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samplesheet = "$projectDir/data/samplesheet_dev.tsv"
params.chrM_fasta = "/references/Homo_sapiens_assembly38.chrM.fasta"
params.chrM_shift_fasta = "/references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta"
params.outdir = "$projectDir/results"
params.shiftback_chain = "/references/ShiftBack.chain"
params.blacklist_sites = "/references/blacklist_sites.hg38.chrM.bed"
params.gnomad = "/references/gnomad.genomes.v3.1.sites.chrM.vcf.bgz"
params.varnote_annoc = "/references/varnote.annoc"


include { MITO_SNVS_INDELS } from './workflows/mito_snvs_indels'

workflow {
    main:

    Channel
        .fromPath(params.samplesheet)
        | splitCsv( header: true , sep: "\t")
        | map { row ->
            meta = row.subMap('sample_id', 'replicate', 'family', 'source' )
            [meta, file(row.bam_file, checkIfExists: true)]
        }
        | set { samplesheet }


    MITO_SNVS_INDELS (
        samplesheet,
        params.chrM_fasta,
        params.chrM_shift_fasta,
        params.shiftback_chain,
        params.blacklist_sites,
        params.gnomad,
        params.varnote_annoc
    )
}