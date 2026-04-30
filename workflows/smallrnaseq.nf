/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_GENOME                   } from '../subworkflows/local/prepare_genome/main'
include { CAT_FASTQ                        } from '../modules/nf-core/cat/fastq/main'
include { FASTQ_TRIM_FASTP_FASTQC          } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { FASTQC as FASTQC_MERGED          } from '../modules/nf-core/fastqc/main'
include { BOWTIE2_BUILD                    } from '../modules/nf-core/bowtie2/build/main'
include { FASTQ_ALIGN_BOWTIE2              } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { BOWTIE_BUILD                     } from '../modules/nf-core/bowtie/build/main'
include { FASTQ_ALIGN_BOWTIE               } from '../subworkflows/local/fastq_align_bowtie/main'
include { SHORTSTACK                       } from '../modules/local/shortstack/main'
include { SHORTSTACK as SHORTSTACK_COMPARE } from '../modules/local/shortstack/main'
include { DESEQ2_QC                        } from '../modules/local/deseq2_qc/main'
include { SHORTSTACK_SUM_COUNTS            } from '../modules/local/shortstack_sum_counts/main'
include { MULTIQC                          } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText           } from '../subworkflows/local/utils_nfcore_smallrnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SMALLRNASEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    /*
    ================================================================================
                                    Prepare reference files
    ================================================================================
    */

    ch_fasta_raw      = params.fasta ? channel.fromPath(params.fasta).map { fasta -> [ [id:file(fasta).baseName], fasta ] }.collect() : channel.empty()
    ch_annotation_raw = params.gtf   ? channel.fromPath(params.gtf).map   { gtf -> [ [id:file(gtf).baseName], gtf ] }.collect() :
                        params.gff   ? channel.fromPath(params.gff).map   { gff -> [ [id:file(gff).baseName], gff ] }.collect() : channel.empty()

    PREPARE_GENOME ( ch_fasta_raw, ch_annotation_raw )

    ch_fasta      = PREPARE_GENOME.out.fasta
    ch_fai        = PREPARE_GENOME.out.fai
    ch_gzi        = PREPARE_GENOME.out.gzi
    ch_annotation = PREPARE_GENOME.out.annotation

    /*
    ================================================================================
                                    FASTQ concatenation
    ================================================================================
    */

    CAT_FASTQ (
        ch_samplesheet.map { meta, fastqs -> [ meta, fastqs ] }
    )

    /*
    ================================================================================
                                    QC and trimming
    ================================================================================
    */

    if (!params.skip_fastp & !params.skip_fastqc) {
        FASTQ_TRIM_FASTP_FASTQC (
            CAT_FASTQ.out.reads.map { meta, reads -> [ meta, reads, [] ] },
            params.save_trimmed_fail,
            params.discard_trimmed_pass,
            params.save_merged,
            params.skip_fastp,
            params.skip_fastqc
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_raw_zip.map { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.trim_json.map { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_FASTP_FASTQC.out.fastqc_trim_zip.map { it[1] })

        ch_trim_reads_merged = FASTQ_TRIM_FASTP_FASTQC.out.trim_reads_merged
        FASTQC_MERGED ( ch_trim_reads_merged )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_MERGED.out.zip.map { it[1] })

    } else {
        ch_trim_reads_merged = CAT_FASTQ.out.reads
    }

    /*
    ================================================================================
                                    Alignment
    ================================================================================
    */

    //
    // MODULE: Bowtie2 build
    //
    BOWTIE2_BUILD ( ch_fasta.join(ch_fai) )

    ch_bowtie2_input = ch_trim_reads_merged
        .map { meta, reads -> [ meta + [ single_end: true ], reads ] }
    FASTQ_ALIGN_BOWTIE2 (
        ch_bowtie2_input,
        BOWTIE2_BUILD.out.index.first(),
        params.save_unaligned,
        false,
        ch_fasta.join(ch_fai).first()
    )
    ch_bowtie2_bam = FASTQ_ALIGN_BOWTIE2.out.bam
        .map { meta, bam -> [ meta + [ aligner: "Bowtie2" ], bam ] }
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.log_out.map { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.stats.map { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.flagstat.map { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.idxstats.map { it[1] })

    //
    // MODULE: Bowtie build
    // Mainly a sanity check to ensure that we get similar results to Bowtie2
    //
    ch_bowtie_bam = channel.empty()
    if (!params.skip_bowtie) {
        BOWTIE_BUILD ( ch_fasta )

        ch_bowtie_input = ch_trim_reads_merged
            .map { meta, reads -> [ meta + [ single_end: true ], reads ] }
        FASTQ_ALIGN_BOWTIE (
            ch_bowtie_input,
            BOWTIE_BUILD.out.index.first(),
            params.save_unaligned,
            ch_fasta.join(ch_fai).first()
        )
        ch_bowtie_bam = FASTQ_ALIGN_BOWTIE.out.bam
            .map { meta, bam -> [ meta + [ aligner: "Bowtie" ], bam ] }
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.log_out.map { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.stats.map { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.flagstat.map { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.idxstats.map { it[1] })
    }

    /*
    ================================================================================
                                    Discover small RNA loci
    ================================================================================
    */

    // Combine results
    ch_bam = ch_bowtie2_bam
        .mix(ch_bowtie_bam)

    SHORTSTACK (
        ch_bam,
        ch_fasta.map { meta, fasta -> fasta }.first()
    )

    ch_shortstack_input = ch_trim_reads_merged
        .map { meta, reads -> reads }
        .collect()
        .map { bams -> [ [id: 'ShortStack.combined'], bams ] }

    // Entire dataset/s doesn't work on HPC... might be something to do with the program.
    SHORTSTACK_COMPARE (
        ch_shortstack_input,
        ch_fasta.map { meta, fasta -> fasta }.first()
    )

    /*
    ================================================================================
                                    Post-alignment QC
    ================================================================================
    */

    ch_biotype_input = SHORTSTACK_COMPARE.out.bam
        .transpose()
        .map { meta, bam -> [ [id: bam.baseName], bam ] }
        .combine(
            SHORTSTACK_COMPARE.out.counts.join(SHORTSTACK_COMPARE.out.gff3)
            .map { meta, counts, gff -> [ counts, gff ] }
        )
        .map { meta, bam, counts, gff -> [ meta, gff, counts ] }

    SHORTSTACK_SUM_COUNTS ( ch_biotype_input )
    ch_multiqc_files = ch_multiqc_files.mix(SHORTSTACK_SUM_COUNTS.out.counts_txt)

    DESEQ2_QC ( SHORTSTACK_COMPARE.out.counts )
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.pca_txt)
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.rle_txt)
    ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.dists_txt)

    //  
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'smallrnaseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'smallrnaseq'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
