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
                | SHORTEN_SEQ_IDS_IF_REQ
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

process SHORTEN_SEQ_IDS_IF_REQ {
    tag "${hap_name}"
    label "process_single"
    
    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("*.renamed.ids.fasta")
    
    script:
        """
        cat $fasta_file | grep -o '^>[^ ]*' | sed 's/>//1' > input_file_ids.txt
        shorten_fasta_ids_ba0fcb9.py input_file_ids.txt > short_ids.tsv

        fasta_file_bash_var="$fasta_file"
        output_file="\${fasta_file_bash_var%%.*}.renamed.ids.fasta"

        if [[ "\$(cat short_ids.tsv)" =~ "IDs have acceptable length and character" ]];
        then
            echo "IDs have acceptable length and character"
            cat "$fasta_file" > "\$output_file"
            exit 0
        fi

        declare -A substitution_mapping
        while IFS=\$'\\t' read -r fasta_id substitute_id; do
            substitution_mapping["\$fasta_id"]="\$substitute_id"
        done < short_ids.tsv

        while IFS= read -r line; do
            if [[ \$line =~ ^\\>([^[:space:]]+)(.*) ]]; then
                fasta_id="\${BASH_REMATCH[1]}"
                echo ">\${substitution_mapping["\$fasta_id"]}\${BASH_REMATCH[2]}" >> "\$output_file"
            else
                echo "\$line" >> "\$output_file"
            fi
        done < "$fasta_file"
        """
}

process EDTA {
    tag "${hap_name}"
    label "process_high"
    label "process_week_long"
    
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
        --threads ${task.cpus}

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
    tag "${hap_name}"
    if (params.lai.mode != "-qq") {
        label "process_high"
        label "process_week_long"
    } else {
        label "process_single"
    }
    
    container 'quay.io/biocontainers/ltr_retriever:2.9.0--hdfd78af_1'
    publishDir "${params.outdir.main}/lai", mode: 'copy'
    
    input:
        tuple val(hap_name), path(fasta_file), path(pass_list), path(genome_out)
    
    output:
        tuple path('*.LAI.log'), path('*.LAI.out') // reversed file name to avoid conflict 
    
    script:
        """
        LAI ${params.lai.mode} \
        -t ${task.cpus} \
        -genome $fasta_file \
        -intact $pass_list \
        -all $genome_out > "${hap_name}.LAI.log"

        lai_output_file_name="\$(basename $pass_list .pass.list).out.LAI"
        [[ -f "\$lai_output_file_name" ]] && cat "\$lai_output_file_name" > "${hap_name}.LAI.out" || echo "LAI OUTPUT IS EMPTY" > "${hap_name}.LAI.out"
        """
}