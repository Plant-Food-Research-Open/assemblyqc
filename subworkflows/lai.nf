nextflow.enable.dsl=2

workflow LAI {
    take:
        tuple_of_hap_genome_pass_out
    
    main:
        if (!params.lai.skip) {
            if (params.lai.pass_list == null || params.lai.out_file == null) {
                tuple_of_hap_genome_pass_out
                | map {
                    return [it[0], it[1]] // [tag, genome fasta file path]
                }
                | EDTA
                | map {
                    return [it[0], it[1], it[3], it[4]] // [tag, genome fasta file path, pass list path, out file path]
                }
                | RUN_LAI
                | collect
                | set { ch_list_of_lai_outputs }
            } else {
                tuple_of_hap_genome_pass_out
                | RUN_LAI
                | collect
                | set { ch_list_of_lai_outputs }
            }
        } else {
            ch_list_of_lai_outputs = Channel.of([])
        }
    
    emit:
        list_of_outputs = ch_list_of_lai_outputs
}

process EDTA {
    label 'uses_high_cpu_mem'
    label 'takes_six_days'
    tag "${hap_name}"
    container 'quay.io/biocontainers/edta:2.1.0--hdfd78af_1'

    publishDir "${params.outdir.main}/edta", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path('*.EDTA.fasta'), path('*.EDTA.TEanno.gff3'), path('*.EDTA.pass.list'), path('*.EDTA.out'), path('*.EDTA.TElib.fa')
    
    script:
        """
        EDTA.pl \
        --genome $fasta_file \
        --step all \
        --sensitive ${params.lai.edta.is_sensitive} \
        --anno 1 \
        --force 1 \
        --species "others" \
        --threads ${task.cpus * params.ht_factor}

        fasta_file_var="$fasta_file"
        fasta_file_base_name="\${fasta_file_var%.*}"
        edta_mod_str=".mod"
        fasta_file_mod="\${fasta_file_var}\${edta_mod_str}"
        
        ln -s "\$fasta_file_mod" "\${fasta_file_base_name}.EDTA.fasta"
        
        [[ -f "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list" ]] \
        && echo "EDTA pass list detected" \
        || echo "EDTA PASS LIST IS EMPTY" > "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list"
        
        ln -s "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list" "\${fasta_file_base_name}.EDTA.pass.list"
        
        ln -s "./\${fasta_file_mod}.EDTA.anno/\${fasta_file_mod}.out" "\${fasta_file_base_name}.EDTA.out"
        """
}

process RUN_LAI {
    label 'uses_high_cpu_mem'
    if (params.lai.mode != "-qq") {
        label 'takes_six_days'
    }
    tag "${hap_name}"
    container 'quay.io/biocontainers/ltr_retriever:2.9.0--hdfd78af_1'

    publishDir "${params.outdir.main}/lai", mode: 'copy'
    
    input:
        tuple val(hap_name), path(genome_fasta), path(pass_list), path(genome_out)
    
    output:
        tuple path('*.LAI.log'), path('*.LAI.out') // reversed file name to avoid conflict 
    
    script:
        """
        LAI ${params.lai.mode} \
        -t ${task.cpus * params.ht_factor} \
        -genome $genome_fasta \
        -intact $pass_list \
        -all $genome_out > "${hap_name}.LAI.log"

        lai_output_file_name="\$(basename $pass_list .pass.list).out.LAI"
        [[ -f "\$lai_output_file_name" ]] && cat "\$lai_output_file_name" > "${hap_name}.LAI.out" || echo "LAI OUTPUT IS EMPTY" > "${hap_name}.LAI.out"
        """
}