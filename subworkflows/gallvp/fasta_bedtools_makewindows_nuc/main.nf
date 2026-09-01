include { SAMTOOLS_FAIDX            } from '../../../modules/gallvp/samtools/faidx/main'
include { BEDTOOLS_MAKEWINDOWS      } from '../../../modules/gallvp/bedtools/makewindows/main'
include { BEDTOOLS_NUC              } from '../../../modules/gallvp/bedtools/nuc/main'

workflow FASTA_BEDTOOLS_MAKEWINDOWS_NUC {

    take:
    ch_fasta                        // channel: [ val(meta), fasta ]

    main:

    // MODULES: SAMTOOLS_FAIDX
    SAMTOOLS_FAIDX (
        ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] },
        true // get_sizes
    )

    ch_sizes                        = SAMTOOLS_FAIDX.out.sizes

    // collectFile: Regions BED
    ch_regions_bed                  = ch_sizes
                                    | map { meta, sizes ->
                                        [
                                            "${meta.id}.regions.bed",
                                            sizes.readLines().collect { row ->
                                                def cols = row.tokenize('\t')
                                                def name = cols[0]
                                                def size = cols[1]
                                                def start = 0
                                                def end = size // BED is 0-indexed end-exclusive

                                                "${name}\t${start}\t${end}"
                                            }.join('\n')
                                        ]
                                    }
                                    | collectFile(newLine: true)
                                    | map { file ->
                                        [
                                            [ id: file.baseName.replace('.regions', '') ], // meta2: [ id: meta.id ]
                                            file
                                        ]
                                    }

    // MODULE: BEDTOOLS_MAKEWINDOWS
    BEDTOOLS_MAKEWINDOWS (
        ch_regions_bed
    )

    ch_intervals_bed                = BEDTOOLS_MAKEWINDOWS.out.bed

    // MODULE: BEDTOOLS_NUC
    ch_bedtools_nuc_inputs          = ch_fasta
                                    | map { meta, fasta -> [ [ id: meta.id ], fasta ] }
                                    | join(
                                        ch_intervals_bed
                                    )
    BEDTOOLS_NUC(
        ch_bedtools_nuc_inputs
    )

    emit:
    fai                             = SAMTOOLS_FAIDX.out.fai    // channel: [ val(meta), fai ]
    bed                             = ch_intervals_bed          // channel: [ val(meta), bed ]
    nuc                             = BEDTOOLS_NUC.out.bed      // channel: [ val(meta2), bed ]
}
