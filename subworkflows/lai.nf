nextflow.enable.dsl=2

workflow LAI {
    take:
        tuple_of_hap_file
    
    main: 
        EDTA(tuple_of_hap_file)
        | RUN_LAI
        | collect
        | set { ch_list_of_lai_logs }
    
    emit:
        list_of_lai_logs = ch_list_of_lai_logs
}

process EDTA {
    label 'uses_high_cpu_mem'
    label 'takes_days'
    tag "${hap_name}"
    container 'quay.io/biocontainers/edta:2.1.0--hdfd78af_1'

    publishDir "${params.outdir.main}/edta", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        val hap_name, emit: hap_name
        path '*.EDTA.fasta', emit: genome_fasta
        path '*.EDTA.TEanno.gff3', emit: te_anno_gff3
        path '*.EDTA.pass.list', emit: pass_list
        path '*.EDTA.out', emit: genome_out
    
    script:
        """
        EDTA.pl \
        --genome $fasta_file \
        --step all \
        --sensitive ${params.lai.edta.is_sensitive} \
        --anno 1 \
        --force 1 \
        --species "others" \
        --threads ${task.cpus}

        fasta_file_var="$fasta_file"
        fasta_file_base_name="\${fasta_file_var%.*}"
        edta_mod_str=\$(cat .command.log | grep "Trying to reformat seq IDs..." > \\dev\\null && echo ".mod" || echo "")
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
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    container 'quay.io/biocontainers/ltr_retriever:2.9.0--hdfd78af_1'
    
    input:
        val hap_name
        path genome_fasta
        path te_anno_gff3
        path pass_list
        path genome_out
    
    output:
        path '*.LAI.log'
    
    script:
        """
        LAI \
        -genome $genome_fasta \
        -intact $pass_list \
        -all $genome_out > "${hap_name}.LAI.log"
        """
}