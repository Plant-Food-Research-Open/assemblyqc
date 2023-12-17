nextflow.enable.dsl=2

include { GZIP_FASTA } from '../../modules/local/gzip_fasta'

workflow SYNTENY {
    take:
        tuple_of_tag_fasta_seq_list
        tuple_of_tag_xref_fasta_seq_list
    
    main:
        if(!params.synteny.skip) {
        
            if(params.synteny.between_target_asm) {
                tuple_of_tag_fasta_seq_list
                | map {
                    [it]
                }
                | collect
                | map {
                    getUniqueWithinCombinations(it)
                }
                | flatten
                | buffer(size:6)
                | set { ch_between_target_asm_combinations }
            } else {
                ch_between_target_asm_combinations = Channel.empty()
            }

            tuple_of_tag_xref_fasta_seq_list
            | map {
                [it[0], it[1]] // [tag, xref fasta file path]
            }
            | GZIP_FASTA
            | join(
                tuple_of_tag_xref_fasta_seq_list
            )
            | map {
                [it[0], it[1], it[3]] // [tag, uncompressed xref fasta file path, seq list]
            }
            | set { ch_tuple_tag_xref_uncompressed_fasta_seq_list }
            
            ch_between_target_asm_combinations
            .mix(
                tuple_of_tag_fasta_seq_list
                | combine(
                    ch_tuple_tag_xref_uncompressed_fasta_seq_list
                )
            )
            .tap { ch_full_tap_from_all_combinations }
            .map {
                ["${it[0]}.on.${it[3]}", it[2], it[5]] // [target.on.reference, target_seq_list, ref_seq_list]    
            }
            .set { ch_seq_lists }
            
            
            ch_full_tap_from_all_combinations
            | FILTER_SORT_FASTA_AND_VALIDATE_SEQ_LISTS
            | (MUMMER & GET_FASTA_LEN)
            
            MUMMER
            .out
            .tag_delta_file
            | DNADIFF
            | CIRCOS_BUNDLE_LINKS
            | ADD_COLOUR_TO_BUNDLE_LINKS
            | join(ch_seq_lists)
            | RELABEL_BUNDLE_LINKS
            | SPLIT_BUNDLE_FILE_BY_TARGET_SEQS
            | map {
                flattenSplitBundles(it)
            }
            | flatten
            | buffer(size:3)
            | set { ch_circos_split_bundle_links }

            GET_FASTA_LEN
            .out
            .tag_len_files
            | join(ch_seq_lists)
            | RELABEL_FASTA_LEN
            | cross(
                ch_circos_split_bundle_links
            )
            | map {
                [it[0][0], it[1][1], it[1][2], it[0][1], it[0][2]] // [target.on.reference, seq_tag, split_bundle_file, target_seq_len, ref_seq_len]
            }
            | GENERATE_KARYOTYPE
            | join(
                ch_circos_split_bundle_links
                | map {
                    ["${it[0]}.${it[1]}", it[2]] // [target.on.reference.seq_tag, split_bundle_file]
                }
            )
            | CIRCOS

            CIRCOS
            .out
            .png_file
            | collect
            | set{ ch_list_of_circos_plots }
        }
        else {
            ch_list_of_circos_plots = Channel.of([])
        }
    
    emit:
        list_of_circos_plots = ch_list_of_circos_plots
}

def getUniqueWithinCombinations(inputArray) {
    if (inputArray.size() <= 1) {
        return []
    }

    inputArray.sort { a, b -> a[0].compareTo(b[0]) }

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

def flattenSplitBundles(inputArray) {
    def target_on_ref = inputArray[0]
    def files = inputArray[1]

    if(files in ArrayList) {
        return files.collect { [target_on_ref, extractBundleTag(it), it] }
    } else {
        return [files].collect { [target_on_ref, extractBundleTag(it), it] }
    }
}

def extractBundleTag(filePath) {
    def regex = /.*\.(\w+)\.split\.bundle\.txt/
    def matcher = filePath =~ regex
    if (matcher.matches()) {
        return matcher.group(1)
    } else {
        // This branch should never be executed if all the upstream logic is implemented correctly.
        error "Error: Failed to parse the sequence tag from file name: ${filePath.getName()}"
    }
}

process FILTER_SORT_FASTA_AND_VALIDATE_SEQ_LISTS {
    tag "${target}.on.${reference}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1"

    input:
        tuple val(target), path(target_fasta), path(target_seq_list), val(reference), path(ref_fasta), path(ref_seq_list)
    
    output:        
        tuple val(target), val(reference), path("filtered.ordered.target.fasta"), path("filtered.ordered.ref.fasta"), emit: tags_fasta_files
    
    script:
        """
        validate_seq_lists_1d50376.sh "$target_seq_list" "$ref_seq_list"
        samtools faidx $target_fasta \$(awk '{print \$1}' $target_seq_list) > filtered.ordered.target.fasta
        samtools faidx $ref_fasta \$(awk '{print \$1}' $ref_seq_list) > filtered.ordered.ref.fasta
        """
}

process GET_FASTA_LEN {
    tag "${target}.on.${reference}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1"
    
    input:        
        tuple val(target), val(reference), path(filtered_ordered_target_fasta), path(filtered_ordered_ref_fasta)
    
    output:
        tuple val("${target}.on.${reference}"), path("target.seq.lengths"), path("ref.seq.lengths"), emit: tag_len_files
    
    script:
        """
        samtools faidx $filtered_ordered_target_fasta
        samtools faidx $filtered_ordered_ref_fasta

        cat "${filtered_ordered_target_fasta}.fai" | awk '{print \$1, \$2}' OFS="\\t" > target.seq.lengths
        cat "${filtered_ordered_ref_fasta}.fai" | awk '{print \$1, \$2}' OFS="\\t" > ref.seq.lengths
        """
}

process MUMMER {
    tag "${target}.on.${reference}"
    label "process_high"
    
    container "docker://staphb/mummer:4.0.0"

    input:
        tuple val(target), val(reference), path(target_fasta), path(ref_fasta)
    
    output:
        tuple val("${target}.on.${reference}"), path("*.delta"), emit: tag_delta_file
    
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
    label "process_single"
    label "process_week_long"
    
    container "docker://staphb/mummer:4.0.0"

    input:
        tuple val(target_on_ref), path(dnadiff_file)
    
    output:
        tuple val(target_on_ref), path("*.xcoords"), path("*.report")
    
    script:
        """
        dnadiff \
        -p $target_on_ref \
        -d $dnadiff_file

        if [[ "${params.synteny.many_to_many_align}" = "1" ]];then
            cat "${target_on_ref}.mcoords" > "${target_on_ref}.m.xcoords"
        else
            cat "${target_on_ref}.1coords" > "${target_on_ref}.1.xcoords"
        fi
        """
}

process CIRCOS_BUNDLE_LINKS {
    tag "${target_on_ref}"
    label "process_single"
    
    container "docker://gallvp/circos-tools:v0.23-1_ps"

    input:
        tuple val(target_on_ref), path(coords_file), path(report_file)
    
    output:
        tuple val(target_on_ref), path("*.xcoords.bundle.txt")
    
    script:
        """
        cat $coords_file | awk '{print \$12,\$1,\$2,\$13,\$3,\$4}' OFS="\\t" > "\$(basename $coords_file).links.txt"

        /usr/share/circos/tools/bundlelinks/bin/bundlelinks \
        -links "\$(basename $coords_file).links.txt" \
        -max_gap "${params.synteny.max_gap}" \
        -min_bundle_size "${params.synteny.min_bundle_size}" \
        1>"\$(basename $coords_file).bundle.txt" \
        2>bundlelinks.err
        """
}

process ADD_COLOUR_TO_BUNDLE_LINKS {
    tag "${target_on_ref}"
    label "process_single"
    
    container "docker://gallvp/python3npkgs:v0.4"

    input:
        tuple val(target_on_ref), path(bundle_links)
    
    output:
        tuple val(target_on_ref), path("*.xcoords.bundle.coloured.txt"), emit: coloured_bundle_links
    
    script:
        """
        if [[ "${params.synteny.color_by_contig}" = "1" ]];then
            color_circos_bundles_by_contig_943e0fb.py \
            "${bundle_links}" \
            > "\$(basename $bundle_links .bundle.txt).bundle.coloured.txt"
        else
            add_color_2_circos_bundle_file_943e0fb.pl \
            -i="${bundle_links}" \
            -o="\$(basename $bundle_links .bundle.txt).bundle.coloured.txt"
        fi
        """
}

process RELABEL_BUNDLE_LINKS {
    tag "${target_on_ref}"
    label "process_single"
    
    container "docker://gallvp/python3npkgs:v0.4"
    
    input:
        tuple val(target_on_ref), path(coloured_bundle_links), path(target_seq_list), path(ref_seq_list)
    
    output:
        tuple val(target_on_ref), path("*.xcoords.bundle.coloured.relabeled.txt"), emit: relabeled_coloured_bundle_links
    
    script:
        """
        #!/usr/bin/env python

        import pandas as pd
        import sys
        import os

        output_file_name = ".".join("$coloured_bundle_links".split(".")[0:-1]) + ".relabeled.txt"

        subs_target_seq = pd.read_csv('$target_seq_list', sep='\\t', header=None)
        subs_target_seq_dict = dict(zip(subs_target_seq.iloc[:, 0], subs_target_seq.iloc[:, 1]))

        subs_ref_seq = pd.read_csv('$ref_seq_list', sep='\\t', header=None)
        subs_ref_seq_dict = dict(zip(subs_ref_seq.iloc[:, 0], subs_ref_seq.iloc[:, 1]))
        
        if os.path.getsize('$coloured_bundle_links') == 0:
            with open(output_file_name, 'w') as f:
                f.write('')
            sys.exit(0)
        else:
            df = pd.read_csv('$coloured_bundle_links', sep=' ', header=None)
        
        df.iloc[:, 3] = df.iloc[:, 3].replace(subs_target_seq_dict, regex=False)
        df.iloc[:, 0] = df.iloc[:, 0].replace(subs_ref_seq_dict, regex=False)
        
        df.to_csv(output_file_name, sep=' ', index=False, header=None)
        """
}

process RELABEL_FASTA_LEN {
    tag "${target_on_ref}"
    label "process_single"
    
    container "docker://gallvp/python3npkgs:v0.4"
    
    input:
        tuple val(target_on_ref), path(target_seq_lengths), path(ref_seq_lengths), path(target_seq_list), path(ref_seq_list)
    
    output:
        tuple val(target_on_ref), path("relabeld.target.seq.lengths"), path("relabeld.ref.seq.lengths"), emit: relabeled_seq_lengths
    
    script:
        """
        #!/usr/bin/env python

        import pandas as pd

        subs_target_seq = pd.read_csv('$target_seq_list', sep='\\t', header=None)
        subs_target_seq_dict = dict(zip(subs_target_seq.iloc[:, 0], subs_target_seq.iloc[:, 1]))

        subs_ref_seq = pd.read_csv('$ref_seq_list', sep='\\t', header=None)
        subs_ref_seq_dict = dict(zip(subs_ref_seq.iloc[:, 0], subs_ref_seq.iloc[:, 1]))
        
        df_target_seq_lengths = pd.read_csv('$target_seq_lengths', sep='\\t', header=None)
        df_target_seq_lengths.iloc[:, 0] = df_target_seq_lengths.iloc[:, 0].replace(subs_target_seq_dict, regex=False)
        df_target_seq_lengths.to_csv("relabeld.target.seq.lengths", sep='\\t', index=False, header=None)

        df_ref_seq_lengths = pd.read_csv('$ref_seq_lengths', sep='\\t', header=None)
        df_ref_seq_lengths.iloc[:, 0] = df_ref_seq_lengths.iloc[:, 0].replace(subs_ref_seq_dict, regex=False)
        df_ref_seq_lengths.to_csv("relabeld.ref.seq.lengths", sep='\\t', index=False, header=None)
        """
}

process SPLIT_BUNDLE_FILE_BY_TARGET_SEQS {
    tag "${target_on_ref}"
    label "process_single"

    input:
        tuple val(target_on_ref), path(coloured_bundle_links)
    
    output:
        tuple val(target_on_ref), path("*.split.bundle.txt")
    
    script:
        """
        if [[ "${params.synteny.plot_1_vs_all}" = "1" ]];then
            target_seqs=(\$(awk '{print \$4}' $coloured_bundle_links | sort | uniq))
            
            for i in "\${!target_seqs[@]}"
            do
                target_seq=\${target_seqs[\$i]}
                awk -v seq="\$target_seq" '\$4==seq {print \$0}' $coloured_bundle_links > "${target_on_ref}.\${target_seq}.split.bundle.txt"
            done
        fi

        cat $coloured_bundle_links > "${target_on_ref}.all.split.bundle.txt"
        """
}

process GENERATE_KARYOTYPE {
    tag "${target_on_ref}.${seq_tag}"
    label "process_single"

    input:
        tuple val(target_on_ref), val(seq_tag), path(split_bundle_file), path(target_seq_len), path(ref_seq_len)
    
    output:
        tuple val("${target_on_ref}.${seq_tag}"), path("*.karyotype")
    
    script:
        """
        ref_seqs=(\$(awk '{print \$1}' $split_bundle_file | sort | uniq))

        if [ \${#ref_seqs[@]} -eq 0 ]; then
            touch "${target_on_ref}.${seq_tag}.karyotype"
            exit 0
        fi

        tmp_file=\$(mktemp)
        printf '%s\\n' "\${ref_seqs[@]}" > "\$tmp_file"

        if [[ $seq_tag = "all" ]];then
            cat $target_seq_len > filtered.target.seq.len
        else
            grep -w "$seq_tag" $target_seq_len > filtered.target.seq.len
        fi
        cat filtered.target.seq.len | awk '{print \$1,\$2,"grey"}' OFS="\\t" > colored.filtered.target.seq.len

        grep -w -f "\$tmp_file" $ref_seq_len > filtered.ref.seq.len
        cat filtered.ref.seq.len | awk '{print \$1,\$2,"black"}' OFS="\\t" > colored.filtered.ref.seq.len

        cat colored.filtered.ref.seq.len | sort -k1V > merged.seq.lengths
        cat colored.filtered.target.seq.len | sort -k1Vr >> merged.seq.lengths
        sed -i '/^\$/d' merged.seq.lengths
        
        cat merged.seq.lengths \
        | awk '{print "chr -",\$1,\$1,"0",\$2-1,\$3}' OFS="\\t" \
        > "${target_on_ref}.${seq_tag}.karyotype"

        rm "\$tmp_file"
        """
}

process CIRCOS {
    tag "${target_on_ref_seq}"
    label "process_single"
    
    container "docker://gallvp/circos-tools:v0.23-1_ps" 
    publishDir "${params.outdir}/synteny/${target_on_ref_seq}", mode: 'copy'

    input:
        tuple val(target_on_ref_seq), path(karyotype), path(bundle_file)
    
    output:
        path "*.svg", emit: svg_file
        path "*.png", emit: png_file
        path "bundled.links.tsv", emit: bundled_links_tsv
        path "circos.conf", emit: circos_conf
        path "karyotype.tsv", emit: karyotype_tsv
    
    script:
        """

        links_count=\$(wc -l < "$bundle_file")
        max_links=20000
        if [ "\$links_count" -gt "\$max_links" ]; then
            echo "Link count exceeded \$max_links for ${bundle_file}."
            echo "Try to shrink the number of links by increasing the max_gap and min_bundle_size options in the config file."
            exit 1
        fi

        cat $karyotype > "karyotype.tsv"
        cat $bundle_file | awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' OFS="\\t" > bundled.links.tsv

        num_sequences=\$(cat $karyotype | wc -l)
        if (( \$num_sequences <= 10 )); then
            label_font_size=40
        elif (( \$num_sequences <= 30 )); then
            label_font_size=30
        else
            label_font_size=15
        fi

        if (( \$num_sequences <= 10 )); then
            ticks_config="<ticks>
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
            </ticks>"
            
            label_offset=" + 120p"
        else
            ticks_config=""
            
            label_offset=" + 25p"
        fi

        cat <<- EOF > circos.conf
        # circos.conf
        karyotype = $karyotype

        <ideogram>
            <spacing>
                default             = 0.005r
            </spacing>

            radius                  = 0.8r
            thickness               = 25p
            fill                    = yes
            stroke_thickness        = 0

            show_label              = yes
            label_font              = default
            label_radius            = dims(ideogram,radius_outer)\$label_offset
            label_size              = \$label_font_size
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
        
        \$ticks_config
        
        <image>
            <<include /usr/share/circos/etc/image.conf>>
        </image>
        <<include /usr/share/circos/etc/colors_fonts_patterns.conf>>
        <<include /usr/share/circos/etc/housekeeping.conf>>
EOF

        if [ ! -s $karyotype ]; then
            touch "${target_on_ref_seq}.svg"
            touch "${target_on_ref_seq}.png"
            exit 0
        fi

        circos

        mv circos.svg "${target_on_ref_seq}.svg"
        mv circos.png "${target_on_ref_seq}.png"
        """
}