nextflow.enable.dsl=2

include { MAKE_AGP_FROM_FASTA       } from '../modules/make_agp_from_fasta.nf'
include { AGP2_ASSEMBLY             } from '../modules/agp2_assembly.nf'
include { ASSEMBLY2_BEDPE           } from '../modules/assembly2_bedpe.nf'
include { MATLOCK_BAM2_JUICER       } from '../modules/matlock_bam2_juicer.nf'
include { RUN_ASSEMBLY_VISUALIZER   } from '../modules/run_assembly_visualizer.nf'

workflow CREATE_HIC_FILE {
    take:
        assembly_fasta
        bam_file
        hap_tag

    main:
        ch_agp_file                     = MAKE_AGP_FROM_FASTA(assembly_fasta)
        ch_assembly_file                = AGP2_ASSEMBLY(ch_agp_file)
        ch_txt_file                     = MATLOCK_BAM2_JUICER(bam_file)
        ch_hic_file                     = RUN_ASSEMBLY_VISUALIZER(ch_assembly_file, ch_txt_file, hap_tag)
        ch_bedpe_file                   = ASSEMBLY2_BEDPE(ch_assembly_file, hap_tag)

    emit:
        hic_file                        = ch_hic_file
}