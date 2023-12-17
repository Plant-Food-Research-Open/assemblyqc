nextflow.enable.dsl=2

include { MAKE_AGP_FROM_FASTA       } from '../../modules/local/make_agp_from_fasta.nf'
include { AGP2_ASSEMBLY             } from '../../modules/local/agp2_assembly.nf'
include { ASSEMBLY2_BEDPE           } from '../../modules/local/assembly2_bedpe.nf'
include { MATLOCK_BAM2_JUICER       } from '../../modules/local/matlock_bam2_juicer.nf'
include { JUICER_SORT               } from '../../modules/local/juicer_sort.nf'
include { RUN_ASSEMBLY_VISUALIZER   } from '../../modules/local/run_assembly_visualizer.nf'

workflow CREATE_HIC_FILE {
    take:
        create_hic_file_inputs // [sample_id.on.tag, assembly_fasta, alignment_bam]

    main:
        create_hic_file_inputs
        | map {
            [it[0], it[1]] // [sample_id.on.tag, assembly_fasta]
        }
        | MAKE_AGP_FROM_FASTA
        | AGP2_ASSEMBLY
        | ASSEMBLY2_BEDPE

        create_hic_file_inputs
        | map {
            [it[0], it[2]] // [sample_id.on.tag, alignment_bam]
        }
        | MATLOCK_BAM2_JUICER
        | JUICER_SORT

        AGP2_ASSEMBLY.out.agp_assembly_file
        | join(JUICER_SORT.out.sorted_links_txt_file) // [sample_id.on.tag, agp_assembly_file, sorted_links_txt_file]
        | RUN_ASSEMBLY_VISUALIZER

    emit:
        hic_file = RUN_ASSEMBLY_VISUALIZER.out.hic_file // [sample_id_on_tag, hic_file]
}