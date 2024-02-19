include { MAKE_AGP_FROM_FASTA       } from '../../modules/local/make_agp_from_fasta.nf'
include { AGP2_ASSEMBLY             } from '../../modules/local/agp2_assembly.nf'
include { ASSEMBLY2_BEDPE           } from '../../modules/local/assembly2_bedpe.nf'
include { MATLOCK_BAM2_JUICER       } from '../../modules/local/matlock_bam2_juicer.nf'
include { JUICER_SORT               } from '../../modules/local/juicer_sort.nf'
include { RUN_ASSEMBLY_VISUALIZER   } from '../../modules/local/run_assembly_visualizer.nf'

workflow CREATE_HIC_FILE {
    take:
    tuple_of_tag_fa_bam

    main:
    // MODULE: MAKE_AGP_FROM_FASTA | AGP2_ASSEMBLY | ASSEMBLY2_BEDPE
    MAKE_AGP_FROM_FASTA ( tuple_of_tag_fa_bam.map { tag, fa, bam -> [ tag, fa ] } )
    | AGP2_ASSEMBLY
    | ASSEMBLY2_BEDPE

    // MODULE: MATLOCK_BAM2_JUICER | JUICER_SORT
    MATLOCK_BAM2_JUICER ( tuple_of_tag_fa_bam.map { tag, fa, bam -> [ tag, bam ] } )
    | JUICER_SORT

    // MODULE: RUN_ASSEMBLY_VISUALIZER
    RUN_ASSEMBLY_VISUALIZER ( AGP2_ASSEMBLY.out.assembly.join(JUICER_SORT.out.links) )

    emit:
    hic = RUN_ASSEMBLY_VISUALIZER.out.hic
}
