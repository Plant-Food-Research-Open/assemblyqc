include { BWA_INDEX         } from '../../../modules/pfr/bwa/index/main'
include { BWA_MEM           } from '../../../modules/pfr/bwa/mem/main'
include { SAMBLASTER        } from '../../../modules/pfr/samblaster/main'

workflow FASTQ_BWA_MEM_SAMBLASTER {

    take:
    ch_fastq                // channel: [ val(meta), [ fq ] ]
    ch_fasta_bwa_index      // channel: [ fasta, index ]; fast | index

    main:
    ch_versions             = Channel.empty()

    ch_has_index            = ch_fasta_bwa_index
                            | branch { fasta, index ->
                                yes: index
                                no: !index
                            }

    // MODULE: BWA_INDEX
    BWA_INDEX ( ch_has_index.no.map { fasta, index -> [ [], fasta ] } )

    ch_bwa_index            = BWA_INDEX.out.index
                            | mix(
                                ch_has_index.yes
                                | map { fasta, index ->
                                    [ [], index ]
                                }
                            )
                            | map { dummy, index -> index }

    ch_versions             = ch_versions.mix(BWA_INDEX.out.versions.first())

    // MODULE: BWA_MEM
    ch_mem_inputs           = ch_fastq
                            | combine(
                                ch_bwa_index
                            )
    def sort_bam            = false
    BWA_MEM(
        ch_mem_inputs.map { meta, fq, index -> [ meta, fq ] },
        ch_mem_inputs.map { meta, fq, index -> [ [], index ] },
        sort_bam
    )

    ch_mem_bam              = BWA_MEM.out.bam
    ch_versions             = ch_versions.mix(BWA_MEM.out.versions.first())

    // MODULE: SAMBLASTER
    SAMBLASTER ( ch_mem_bam )

    ch_blasted_bam          = SAMBLASTER.out.bam
    ch_versions             = ch_versions.mix(SAMBLASTER.out.versions.first())

    emit:
    bam                     = SAMBLASTER.out.bam    // channel: [ val(meta), bam ]
    versions                = ch_versions           // channel: [ versions.yml ]
}
