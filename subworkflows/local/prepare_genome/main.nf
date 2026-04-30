//
// PREPARE GENOME
//

include { GUNZIP as GUNZIP_FASTA      } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_ANNOTATION } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX              } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {

    take:
    fasta      // channel: [meta, fasta]
    annotation // channel: [meta, annotation]

    main:
    //
    // Unzip FASTA if required
    //
    ch_fasta_branched = fasta.branch {
        meta, file ->
            gz: file.extension == 'gz'
            rest: true
    }

    GUNZIP_FASTA ( ch_fasta_branched.gz )

    ch_fasta_unzipped = GUNZIP_FASTA.out.gunzip.mix(ch_fasta_branched.rest)

    //
    // Index FASTA
    //
    SAMTOOLS_FAIDX ( ch_fasta_unzipped.map { meta, file -> [ meta, file, [] ] }, false )

    //
    // Unzip annotation if required
    //
    ch_annotation_branched = annotation.branch {
        meta, file ->
            gz: file.extension == 'gz'
            rest: true
    }

    GUNZIP_ANNOTATION ( ch_annotation_branched.gz )

    ch_annotation_unzipped = GUNZIP_ANNOTATION.out.gunzip.mix(ch_annotation_branched.rest)

    emit:
    fasta      = ch_fasta_unzipped      // channel: [meta, fasta]
    fai        = SAMTOOLS_FAIDX.out.fai // channel: [meta, fai]
    gzi        = SAMTOOLS_FAIDX.out.gzi // channel: [meta, gzi]
    annotation = ch_annotation_unzipped // channel: [meta, annotation]
}
