

include { PROCESS_INPUT_BAM                 } from '../subworkflows/local/process_input_bam'
include { ALIGN_TO_GENOME as ALIGN_TO_REF   } from '../subworkflows/local/align_to_genome'
include { ALIGN_TO_GENOME as ALIGN_TO_SHIFT } from '../subworkflows/local/align_to_genome'
include { VARIANT_CALLING                   } from '../subworkflows/local/variant_calling'
include { ANNOTATION                        } from '../subworkflows/local/annotation'
include { REPORT                            } from '../subworkflows/local/report'


workflow MITO_SNVS_INDELS {

    take:
        ch_samplesheet          // channel: [ val(meta), path(bam) ]
        val_reference_fasta     // string:  path to reference chrM genome
        val_shift_fasta         // string:  path to (shift) reference chrM genome
        val_shiftback_chain     // string:  path to shiftback chain file
        val_blacklist_sites     // string:  path to variant's black list
        val_gnomad              // string: path to gnomAD VCF database in docker image
        val_annoc               // string: path to VarNote annoc file in docker image

    main:
        PROCESS_INPUT_BAM ( ch_samplesheet )

        ALIGN_TO_REF (
            PROCESS_INPUT_BAM.out.mito_fastq,
            PROCESS_INPUT_BAM.out.revert_bam,
            val_reference_fasta
        )

        ALIGN_TO_SHIFT (
            PROCESS_INPUT_BAM.out.mito_fastq,
            PROCESS_INPUT_BAM.out.revert_bam,
            val_shift_fasta
        )

        VARIANT_CALLING (
            ALIGN_TO_REF.out.sorted,
            ALIGN_TO_SHIFT.out.sorted,
            val_reference_fasta,
            val_shift_fasta,
            val_shiftback_chain,
            val_blacklist_sites
        )

        ANNOTATION (
            VARIANT_CALLING.out.normalized,
            val_gnomad,
            val_annoc
        )

        REPORT (
            ANNOTATION.out.gnomad_annotated,
            ANNOTATION.out.haplocheck,
            ANNOTATION.out.haplogrep,
            ALIGN_TO_REF.out.depth
        )


}