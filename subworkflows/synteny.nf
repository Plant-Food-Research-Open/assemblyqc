nextflow.enable.dsl=2

workflow SYNTENY {
    take:
        tuple_of_hap_genome_seq_list
        tuple_of_tag_file
    
    main:
        if (!params.synteny.skip) {
        
            // Combinations
            tuple_of_hap_genome_seq_list
            | map {
                [it]
            }
            | collect
            | map {
                getUniqueWithinCombinations(it)
            }
            | flatten
            | buffer(size:6)
            | set { ch_within_combinations }

            tuple_of_hap_genome_seq_list
            | combine(
                tuple_of_tag_file
            )
            | set { ch_with_genome_combinations }

            ch_within_combinations
            .mix(ch_with_genome_combinations)
            | FILTER_SORT_FASTA
            | PREFIX_FASTA_IDS
            
            MUMMER(
                PREFIX_FASTA_IDS.out.tags_fasta_files
            )
            | DNADIFF
            | CIRCOS_BUNDLE_LINKS
            | set { ch_circos_bundle_links }
            
            GENERATE_KARYOTYPE(
                PREFIX_FASTA_IDS.out.tags_len_files
            )
            | map {
                appendTags(it[0], it[1].splitText(by: 2))
            }
            | flatten
            | buffer(size:2)
            | set { ch_karyotype }

            ch_circos_bundle_links
            .cross(
                ch_karyotype
            )
            | map {
                [it[0][0], it[0][1], it[0][2], it[1][1]] // [tag.on.tag2, links, bundles, karyotype]
            }
            | CIRCOS

            ch_list_of_synteny_plots = Channel.of([])
        }
        else {
            ch_list_of_synteny_plots = Channel.of([])
        }
    
    emit:
        list_of_synteny_plots = ch_list_of_synteny_plots
}

def getUniqueWithinCombinations(inputArray) {
    if (inputArray.size() <= 1) {
        return []
    }

    def outputList = []

    for (int i = 0; i < inputArray.size() - 1; i++) {
        for (int j = i + 1; j < inputArray.size(); j++) {
            def combination = [
                inputArray[i][0],
                file(inputArray[i][1], checkIfExists: true),
                file(inputArray[i][2], checkIfExists: true),
                inputArray[j][0],
                file(inputArray[j][1], checkIfExists: true),
                file(inputArray[j][2], checkIfExists: true)
            ]
            outputList.add(combination)
        }
    }
    return outputList
}

def appendTags(tag, valuesArray) {
    if (valuesArray.size() <= 1) {
        return []
    }

    def outputList = []

    for (int i = 0; i < valuesArray.size(); i++) {
        outputList.add([tag, valuesArray[i]])
    }
    return outputList
}

process FILTER_SORT_FASTA {
    tag "${target}:${reference}"
    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
        tuple val(target), path(target_fasta), path(target_seq_list), val(reference), path(ref_fasta), path(ref_seq_list)
    
    output:        
        tuple val(target), val(reference), path("filtered.ordered.target.fasta"), path("filtered.ordered.ref.fasta")
    
    script:
        """
        num_target_seq=\$(cat $target_seq_list | wc -l)
        num_ref_seq=\$(cat $ref_seq_list | wc -l)
        min_n=\$((\$num_target_seq>\$num_ref_seq ? \$num_ref_seq : \$num_target_seq))

        head -\$min_n $target_seq_list > target.seq.list
        head -\$min_n $ref_seq_list > reference.seq.list

        samtools faidx $target_fasta \$(cat target.seq.list) > filtered.ordered.target.fasta
        samtools faidx $ref_fasta \$(cat reference.seq.list) > filtered.ordered.ref.fasta
        """
}

process PREFIX_FASTA_IDS {
    tag "${target}:${reference}"
    container "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"
    
    input:        
        tuple val(target), val(reference), path(filtered_ordered_target_fasta), path(filtered_ordered_ref_fasta)
    
    output:
        tuple val(target), val(reference), path("prefixed.filtered.ordered.target.fasta"), path("prefixed.filtered.ordered.ref.fasta"), emit: tags_fasta_files
        tuple val(target), val(reference), path("target.seq.lengths"), path("ref.seq.lengths"), emit: tags_len_files
    
    script:
        """
        cat $filtered_ordered_target_fasta | seqkit replace -p ^ -r "${target}_" > prefixed.filtered.ordered.target.fasta
        cat $filtered_ordered_ref_fasta | seqkit replace -p ^ -r "${reference}_" > prefixed.filtered.ordered.ref.fasta

        awk '/^>/ {if (seqlen){print seqlen}; printf \$1" " ;seqlen=0;next; } { seqlen += length(\$0)}END{print seqlen}' prefixed.filtered.ordered.target.fasta \
        | sed 's/>//1' \
        | awk '{print \$1, \$2}' OFS='\\t' \
        > target.seq.lengths

        awk '/^>/ {if (seqlen){print seqlen}; printf \$1" " ;seqlen=0;next; } { seqlen += length(\$0)}END{print seqlen}' prefixed.filtered.ordered.ref.fasta \
        | sed 's/>//1' \
        | awk '{print \$1, \$2}' OFS='\\t' \
        > ref.seq.lengths
        """
}

process MUMMER {
    tag "${target}.on.${reference}"
    label 'uses_high_cpu_mem'
    label 'uses_64_gb_mem'
    label 'takes_four_hours'
    container "docker://staphb/mummer:4.0.0" 
    publishDir "${params.outdir.main}/synteny/nucmer", mode: 'copy'

    input:
        tuple val(target), val(reference), path(target_fasta), path(ref_fasta)
    
    output:
        tuple val("${target}.on.${reference}"), path("*.delta")
    
    script:
        """
        nucmer \
        --mum \
        -t ${task.cpus} \
        -p "${target}.on.${reference}" \
        $ref_fasta \
        $target_fasta
        """
}

process DNADIFF {
    tag "${target_on_ref}"
    container "docker://staphb/mummer:4.0.0" 
    publishDir "${params.outdir.main}/synteny/dnadiff", mode: 'copy'

    input:
        tuple val(target_on_ref), path(dnadiff_file)
    
    output:
        tuple val(target_on_ref), path("*.1coords"), path("*.report")
    
    script:
        """
        dnadiff \
        -p $target_on_ref \
        -d $dnadiff_file
        """
}

process CIRCOS_BUNDLE_LINKS {
    tag "${target_on_ref}"
    container "docker://gallvp/circos-tools:0.23-1" 
    publishDir "${params.outdir.main}/synteny/bundlelinks", mode: 'copy'

    input:
        tuple val(target_on_ref), path(coords_file), path(report_file)
    
    output:
        tuple val(target_on_ref), path("*.1coords.links.txt"), path("*.1coords.bundle.txt")
    
    script:
        """
        cat $coords_file | awk '{print \$12,\$1,\$2,\$13,\$3,\$4}' OFS="\\t" > "${target_on_ref}.1coords.links.txt"

        /usr/share/circos/tools/bundlelinks/bin/bundlelinks \
        -links "${target_on_ref}.1coords.links.txt" \
        1>"${target_on_ref}.1coords.bundle.txt" \
        2>bundlelinks.err
        """
}

process GENERATE_KARYOTYPE {
    tag "${target}.on.${reference}"

    input:
        tuple val(target), val(reference), path(target_seq_len), path(ref_seq_len)
    
    output:
        tuple val("${target}.on.${reference}"), path("*.karyotype")
    
    script:
        """
        paste -d "\\n" $target_seq_len $ref_seq_len > merged.seq.lengths
        cat merged.seq.lengths | awk '{print "chr -",\$1,\$1,"0",\$2-1,(NR%2==1?"red":"blue")}' OFS="\t" > "${target}.on.${reference}.karyotype"
        """
}

process CIRCOS {
    container "docker://gallvp/circos-tools:0.23-1" 
    publishDir "${params.outdir.main}/synteny", mode: 'copy'

    input:
        tuple val(target_on_ref), path(links_file), path(bundle_file), val(karyotype)
    
    output:
        tuple val(target_on_ref), path("*.svg"), path("*.png")
    
    script:
        """
        cat $bundle_file | awk '{print \$1,\$2,\$3,\$4,\$5,\$6}' OFS="\\t" > bundled.links.tsv
        echo -n "$karyotype" | tee karyotype.tsv
        target_tag=\$(cat karyotype.tsv | awk 'NR==1 {print \$3}')
        ref_tag=\$(cat karyotype.tsv | awk 'NR==2 {print \$3}')

        cat <<- EOF > circos.conf
        # circos.conf
        karyotype = karyotype.tsv

        <ideogram>
            <spacing>
                default             = 0.005r
            </spacing>

            radius                  = 0.90r
            thickness               = 25p
            fill                    = yes
            stroke_thickness        = 0

            show_label              = yes
            label_font              = default 
            label_radius            = dims(ideogram,radius_outer) + 100p
            label_size              = 50
            label_parallel          = yes
        </ideogram>

        <links>
            radius                  = 0.99r
            crest                   = 1
            ribbon                  = yes
            flat                    = yes
            stroke_thickness        = 0
            color                   = grey_a3

            bezier_radius           = 0r
            bezier_radius_purity    = 0.5
            <link>
                file                = bundled.links.tsv
            </link>
        </links>

        show_ticks                  = yes
        show_tick_labels            = yes
        chromosomes_units           = 1000000
        chromosomes_display_default = yes
        <ticks>
        radius                      = dims(ideogram,radius_outer)
        orientation                 = out
        label_multiplier            = 1e-6
        color                       = black
        thickness                   = 5p
        label_offset                = 5p
        <tick>
            spacing                 = 0.5u
            size                    = 10p
            show_label              = yes
            label_size              = 20p
            format                  = %.1f
        </tick>
        <tick>
            spacing                 = 1.0u
            size                    = 15p
            show_label              = yes
            label_size              = 30p
            format                  = %.1f
        </tick>
        </ticks>
        
        <image>
            <<include /usr/share/circos/etc/image.conf>>
        </image>
        <<include /usr/share/circos/etc/colors_fonts_patterns.conf>>
        <<include /usr/share/circos/etc/housekeeping.conf>>
EOF

        circos

        mv circos.svg "\${target_tag}.on.\${ref_tag}.svg"
        mv circos.png "\${target_tag}.on.\${ref_tag}.png"
        """
}

// show_ticks                  = yes
//         show_tick_labels            = yes

//         <ticks>
//             radius                  = 1r
//             color                   = black
//             thickness               = 2p
//             multiplier              = 1e-6
//             format                  = %d

//             <tick>
//                 spacing             = 5u
//                 size                = 10p
//             </tick>

//             <tick>
//                 spacing             = 25u
//                 size                = 15p
//                 show_label          = yes
//                 label_size          = 20p
//                 label_offset        = 10p
//                 format              = %d
//             </tick>
//         </ticks>