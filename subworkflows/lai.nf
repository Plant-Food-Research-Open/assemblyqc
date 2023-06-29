nextflow.enable.dsl=2

workflow LAI {
    take:
        tuple_of_tag_genome_pass_out_mono_seqs
    
    main:
        if (!params.lai.skip) {
            if (params.lai.pass_list.isEmpty() && params.lai.out_file.isEmpty()) {
                tuple_of_tag_genome_pass_out_mono_seqs
                | multiMap {
                    tag_fasta: [it[0], it[1]] // [tag, genome fasta file path]
                    tag_mono_seqs: [it[0], it[4]] // [tag, path for the list of monoploid seqs]
                }
                | set { ch_inputs }
                
                ch_inputs.tag_fasta
                | SHORTEN_SEQ_IDS_IF_REQ
                | multiMap {
                    tag_renamed_fasta: [it[0], it[1]] // [tag, path for the genome fasta file with renamed ids]
                    tag_renamed_ids: [it[0], it[2]] // [tag, path for the list of renamed ids]
                }
                | set { ch_output_shorten_ids }

                ch_output_shorten_ids.tag_renamed_ids
                | join(
                    ch_inputs.tag_mono_seqs
                )
                | map {
                    prepareForMappingMonoSeqs(it)
                }
                | MAP_MONOPLOID_SEQS_TO_NEW_IDS

                ch_output_shorten_ids.tag_renamed_fasta
                | EDTA
                | join(
                    ch_output_shorten_ids.tag_renamed_fasta
                )
                | map {
                    [it[0], it[6], it[3], it[4]] // [tag, path for the genome fasta file with renamed ids, pass list, out file]
                }
                | join(
                    MAP_MONOPLOID_SEQS_TO_NEW_IDS.out.tag_new_ids, remainder: true
                )
                | map {
                    prepareForLAI(it)
                }
                | RUN_LAI
                | collect
                | set { ch_list_of_lai_outputs }
            } else {
                tuple_of_tag_genome_pass_out_mono_seqs
                | map {
                    prepareForLAI(it)
                }
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

def prepareForLAI(inputTuple) {
    if (inputTuple[4] == null) {
        return [inputTuple[0], inputTuple[1], inputTuple[2], inputTuple[3], []]
    } else {
        return inputTuple
    }
}

def prepareForMappingMonoSeqs(inputTuple) {
    if (inputTuple[2] == null) {
        return [inputTuple[0], inputTuple[1], []]
    } else {
        return inputTuple
    }
}

process SHORTEN_SEQ_IDS_IF_REQ {
    tag "${hap_name}"
    label "process_single"

    container "docker://gallvp/python3npkgs:v0.3"
    publishDir "${params.outdir.main}/edta", mode: 'copy', pattern: '*.renamed.ids.tsv'
    
    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("*.renamed.ids.fa"), path("*.renamed.ids.tsv")
    
    script:
        """
        fasta_file_bash_var="$fasta_file"
        output_prefix="\${fasta_file_bash_var%.*}"

        shorten_fasta_ids_c97537f.py "$fasta_file" "\$output_prefix"
        """
}

process EDTA {
    tag "${hap_name}"
    label "process_high"
    label "process_week_long"
    
    container 'quay.io/biocontainers/edta:2.1.0--hdfd78af_1'
    containerOptions "-B $TMPDIR:$TMPDIR"
    publishDir "${params.outdir.main}/edta", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name),
        path('*.EDTA.TEanno.gff3'),
        path('*.EDTA.intact.gff3'),
        path('*.EDTA.pass.list'),
        path('*.EDTA.out'),
        path('*.EDTA.TElib.fa'), emit: edta_outputs
    
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

        fasta_file_mod="${fasta_file}.mod"
        
        [[ -f "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list" ]] \
        && echo "EDTA pass list detected" \
        || echo "EDTA PASS LIST IS EMPTY" > "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list"
        
        ln -s "./\${fasta_file_mod}.EDTA.raw/LTR/\${fasta_file_mod}.pass.list" "\${fasta_file_mod}.EDTA.pass.list"
        
        ln -s "./\${fasta_file_mod}.EDTA.anno/\${fasta_file_mod}.out" "\${fasta_file_mod}.EDTA.out"
        """
}

process MAP_MONOPLOID_SEQS_TO_NEW_IDS {
    tag "${tag_name}"
    label "process_single"

    input:
        tuple val(tag_name), path(renamed_ids_tsv), path(monoploid_seqs)
    
    output:
        tuple val(tag_name), path("*.new.monoploid.seqs.txt"), optional: true, emit: tag_new_ids
    
    script:
        if (!monoploid_seqs.isEmpty())
            """
            cp $monoploid_seqs "${tag_name}.new.monoploid.seqs.txt"
            
            while IFS=\$'\\t' read -r original_id renamed_id; do
                sed -i "s/\$original_id/\$renamed_id/g" "${tag_name}.new.monoploid.seqs.txt"
            done < "$renamed_ids_tsv"
            """
        else
            """
                echo "Empty monoploid seqs file"
            """
}

process RUN_LAI {
    tag "${hap_name}"
    if (params.lai.mode != "-qq") {
        label "process_high"
    } else {
        label "process_single"
    }
    
    container 'quay.io/biocontainers/ltr_retriever:2.9.0--hdfd78af_1'
    publishDir "${params.outdir.main}/lai", mode: 'copy'
    
    input:
        tuple val(hap_name), path(fasta_file), path(pass_list), path(genome_out), path(monoploid_seqs)
    
    output:
        tuple path('*.LAI.log'), path('*.LAI.out') // reversed file name to avoid conflict 
    
    script:
        """
        [[ -z "$monoploid_seqs" ]] && mono_param="" || mono_param="-mono $monoploid_seqs"
        
        LAI ${params.lai.mode} "\$mono_param" \
        -t ${task.cpus} \
        -genome $fasta_file \
        -intact $pass_list \
        -all $genome_out > "${hap_name}.LAI.log"

        lai_output_file_name="\$(basename $pass_list .pass.list).out.LAI"
        [[ -f "\$lai_output_file_name" ]] && cat "\$lai_output_file_name" > "${hap_name}.LAI.out" || echo "LAI OUTPUT IS EMPTY" > "${hap_name}.LAI.out"
        """
}