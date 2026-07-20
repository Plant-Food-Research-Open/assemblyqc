workflow FASTQ_BWA_MEM_SAMBLASTER {

    take:
    ch_fastq
    ch_reference
    val_sort_bam

    main:
    ch_versions             = Channel.empty()

    ch_has_index            = ch_reference
                            | branch { _meta2, _fasta, index ->
                                yes: index
                                no: !index
                            }

    // MODULE: MINIBWA_INDEX
    MINIBWA_INDEX ( ch_has_index.no.map { meta2, fasta, _index -> [ meta2, fasta ] } )

    ch_bwa_index            = MINIBWA_INDEX.out.index
                            | mix(
                                ch_has_index.yes
                                | map { meta2, _fasta, index -> [ meta2, index ] }
                            )


    // MODULE: MINIBWA_MAP
    ch_mem_inputs           = ch_fastq
                            | combine(ch_bwa_index)
                            | map { meta, fq, meta2, index ->
                                [ meta + [ ref_id: meta2.id ], fq, index ]
                            }

    MINIBWA_MAP(
        ch_mem_inputs.map { meta, fq, _index -> [ meta, fq ] },
        ch_mem_inputs.map { _meta, _fq, index -> [ [], index ] },
        [ [], [] ],
        val_sort_bam
    )

    ch_mem_bam              = MINIBWA_MAP.out.aligned

    // MODULE: SAMBLASTER 
    SAMBLASTER ( ch_mem_bam )
    ch_versions             = ch_versions.mix(SAMBLASTER.out.versions.first())

    emit:
    bam                     = SAMBLASTER.out.bam
    versions                = ch_versions    // channel: [ versions.yml ] — old-style only, from SAMBLASTER
}