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
                | map {
                    [it[0], it[1]] // [tag, path for the genome fasta file with renamed ids]
                }
                | set { ch_tag_enamed_ids_fasta }
                
                EDTA(ch_tag_enamed_ids_fasta)
                | join(
                    ch_tag_enamed_ids_fasta
                )
                | map {
                    [it[0], it[6], it[3], it[4]] // [tag, path for the genome fasta file with renamed ids, pass list, out file]
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

    container "docker://gallvp/python3npkgs:v0.2"
    publishDir "${params.outdir.main}/edta", mode: 'copy', pattern: '*.renamed.ids.tsv'
    
    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("*.renamed.ids.fa"), path("*.renamed.ids.tsv")
    
    script:
        """
        fasta_file_bash_var="$fasta_file"
        output_prefix="\${fasta_file_bash_var%.*}"

        status=\$(shorten_fasta_ids_c97537f.py "$fasta_file" "\$output_prefix")

        if [[ "\$status" =~ "IDs have acceptable length and character. No change required." ]];
        then
            echo "\$status" > "\${output_files_prefix}.renamed.ids.tsv"
            ln -s "$fasta_file" "\${output_files_prefix}.renamed.ids.fa"
        fi
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