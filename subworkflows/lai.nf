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

                EDTA.out.edta_outputs
                | join(
                    ch_output_shorten_ids.tag_renamed_ids
                )
                | SAVE_EDTA_OUTPUTS
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
        return [inputTuple[0], inputTuple[1], inputTuple[2], inputTuple[3], [inputTuple[4]]]
    }
}

def prepareForMappingMonoSeqs(inputTuple) {
    if (inputTuple[2] == null) {
        return [inputTuple[0], inputTuple[1], []]
    } else {
        return [inputTuple[0], inputTuple[1], [inputTuple[2]]]
    }
}

process SHORTEN_SEQ_IDS_IF_REQ {
    tag "${tag_name}"
    label "process_single"

    container "docker://gallvp/python3npkgs:v0.4"
    
    input:
        tuple val(tag_name), path(fasta_file)
    
    output:
        tuple val(tag_name), path("*.renamed.ids.fa"), path("*.renamed.ids.tsv")
    
    script:
        """
        fasta_file_bash_var="$fasta_file"
        output_prefix="\${fasta_file_bash_var%.*}"

        shorten_fasta_ids_c97537f.py "$fasta_file" "\$output_prefix"
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
        """
        if [[ -z "$monoploid_seqs" ]]; then
            exit 0
        fi

        renamed_ids_head=\$(head -n 1 "$renamed_ids_tsv")

        if [[ \$renamed_ids_head == "IDs have acceptable length and character. No change required." ]]; then
            cat $monoploid_seqs > "${tag_name}.new.monoploid.seqs.txt"
            exit 0
        fi

        declare -A orig_to_new_ids
        while IFS=\$'\\t' read -r original_id renamed_id; do
            orig_to_new_ids["\$original_id"]="\$renamed_id"
        done < "$renamed_ids_tsv"
        
        while IFS= read -r original_id; do
            echo "\${orig_to_new_ids[\$original_id]}"
        done < "$monoploid_seqs" > "${tag_name}.new.monoploid.seqs.txt"
        """
}

process EDTA {
    tag "${tag_name}"
    label "process_high"
    label "process_week_long"
    
    container 'https://depot.galaxyproject.org/singularity/edta:2.1.0--hdfd78af_1'
    containerOptions "-B $TMPDIR:$TMPDIR"

    input:
        tuple val(tag_name), path(fasta_file)
    
    output:
        tuple val(tag_name),
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

process SAVE_EDTA_OUTPUTS {
    tag "${tag_name}"
    label "process_single"

    publishDir "${params.outdir.main}/edta", mode: 'copy'

    input:
        tuple val(tag_name),
        path(te_anno_gff3),
        path(intact_gff3),
        path(pass_list),
        path(out_file),
        path(te_lib_fa),
        path(renamed_ids_tsv)
    
    output:
        tuple val(tag_name),
        path("${tag_name}.EDTA.TEanno.gff3"),
        path("${tag_name}.EDTA.intact.gff3"),
        path("${tag_name}.renamed.ids.EDTA.pass.list"),
        path("${tag_name}.renamed.ids.EDTA.out"),
        path("${tag_name}.EDTA.TElib.fa"), 
        path("${tag_name}.renamed.ids.tsv"), emit: saved_edta_outputs

    script:
        """
        cat $pass_list > "${tag_name}.renamed.ids.EDTA.pass.list"
        cat $out_file > "${tag_name}.renamed.ids.EDTA.out"
        cat $te_lib_fa > "${tag_name}.EDTA.TElib.fa"
        cat $renamed_ids_tsv > "${tag_name}.renamed.ids.tsv"
        
        renamed_ids_head=\$(head -n 1 "$renamed_ids_tsv")
        
        if [[ \$renamed_ids_head == "IDs have acceptable length and character. No change required." ]]; then
            cat $te_anno_gff3 > "${tag_name}.EDTA.TEanno.gff3"
            cat $intact_gff3 > "${tag_name}.EDTA.intact.gff3"
        else
            reverse_edta_naming_f1b7bce.py "$renamed_ids_tsv" "$te_anno_gff3" "$intact_gff3" "$tag_name"
        fi
        """
}

process RUN_LAI {
    tag "${tag_name}"
    if (params.lai.mode != "-qq") {
        label "process_high"
    } else {
        label "process_single"
    }
    
    container 'https://depot.galaxyproject.org/singularity/ltr_retriever:2.9.0--hdfd78af_1'
    publishDir "${params.outdir.main}/lai", mode: 'copy'
    
    input:
        tuple val(tag_name), path(fasta_file), path(pass_list), path(genome_out), path(monoploid_seqs)
    
    output:
        tuple path('*.LAI.log'), path('*.LAI.out') // reversed file name to avoid conflict 
    
    script:
        """
        if [[ -z "$monoploid_seqs" ]];then
            mono_param=""
            lai_output_file_name="${genome_out}.LAI"
        else
            mono_param="-mono $monoploid_seqs"
            lai_output_file_name="${genome_out}.${monoploid_seqs}.out.LAI"
        fi

        # Remove comments from genome fasta
        sed '/^>/ s/\\s.*\$//' $fasta_file > for.lai.no.comments.fsa
        
        LAI ${params.lai.mode} "\$mono_param" \
        -t ${task.cpus} \
        -genome for.lai.no.comments.fsa \
        -intact $pass_list \
        -all $genome_out > "${tag_name}.LAI.log"

        [[ -f "\$lai_output_file_name" ]] && cat "\$lai_output_file_name" > "${tag_name}.LAI.out" || echo "LAI OUTPUT IS EMPTY" > "${tag_name}.LAI.out"
        """
}