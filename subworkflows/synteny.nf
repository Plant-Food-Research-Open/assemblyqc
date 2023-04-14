nextflow.enable.dsl=2

workflow SYNTENY {
    take:
        tuple_of_hap_genome_seq_list
        tuple_of_tag_file
    
    main:
        if (!params.synteny.skip) {
        
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
            | mix(
                tuple_of_hap_genome_seq_list
                | combine(
                    tuple_of_tag_file
                )
            )
            | map { validateSeqLists(it) }
            | take (1)
            | FILTER_SORT_RELABEL_FASTA
            | MUMMER
            | DNADIFF
            | CIRCOS_BUNDLE_LINKS
            | ADD_COLOUR_TO_BUNDLE_LINKS
            | SPLIT_BUNDLE_FILE_BY_TARGET_SEQS
            | map {
                flattenSplitBundles(it)
            }
            | flatten
            | buffer(size:3)
            | set { ch_circos_split_bundle_links }

            GET_FASTA_LEN(
                FILTER_SORT_RELABEL_FASTA.out.tags_fasta_files
            )
            | cross(
                ch_circos_split_bundle_links
            )
            | map {
                [it[0][0], it[1][1], it[1][2], it[0][1], it[0][2]] // [tag.on.tag, seq_tag, split_bundle_file, target_seq_len, ref_seq_len]
            }
            | GENERATE_KARYOTYPE
            | join(
                ch_circos_split_bundle_links
                | map {
                    ["${it[0]}.${it[1]}", it[2]]
                }
            )
            | take(1)
            | CIRCOS

            CIRCOS
            .out
            .png_file
            .collect()
            .set{ ch_list_of_circos_plots }
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

def validateSeqLists(inputArray) {

    file1 = inputArray[2]
    file2 = inputArray[5]

    def lines1 = file(file1).readLines()
    def lines2 = file(file2).readLines()

    lines1.each { line ->
        def columns = line.split()
        if (columns.size() != 2) {
            throw new Exception("Error: Sequence file ${file1.getName()} does not have exactly two columns.")
        }
    }

    lines2.each { line ->
        def columns = line.split()
        if (columns.size() != 2) {
            throw new Exception("Error: Sequence file ${file2.getName()} does not have exactly two columns.")
        }
    }
    
    def outputLines = lines1 + lines2
    
    def secondColumn = outputLines.collect { it.split()[1] }
    if (secondColumn.size() != secondColumn.unique().size()) {
        throw new Exception("Error: Duplicate sequence labels detected in second column for pair: ${file1.getName()}, ${file2.getName()}")
    }
    
    return inputArray
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
    def tag_on_tag = inputArray[0]
    def files = inputArray[1]
    return files.collect { [tag_on_tag, extractBundleTag(it), it] }
}

def extractBundleTag(filePath) {
   def regex = /.*\.(\w+)\.split\.bundle\.txt/
   def matcher = filePath =~ regex
   if (matcher.matches()) {
      return matcher.group(1)
   } else {
      throw new Exception("Error: Failed to parse the sequence tag from file name: ${filePath.getName()}")
   }
}

process FILTER_SORT_RELABEL_FASTA {
    tag "${target}:${reference}"
    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
        tuple val(target), path(target_fasta), path(target_seq_list), val(reference), path(ref_fasta), path(ref_seq_list)
    
    output:        
        tuple val(target), val(reference), path("filtered.ordered.target.fasta"), path("filtered.ordered.ref.fasta"), emit: tags_fasta_files
    
    script:
        """
        samtools faidx $target_fasta \$(awk '{print \$1}' $target_seq_list) > filtered.ordered.target.fasta
        samtools faidx $ref_fasta \$(awk '{print \$1}' $ref_seq_list) > filtered.ordered.ref.fasta

        while read -r id label; do
            if grep -q "\$id" "filtered.ordered.target.fasta"; then
                sed -i "s/\$id/\$label/g" "filtered.ordered.target.fasta"
            fi
        done < "$target_seq_list"

        while read -r id label; do
            if grep -q "\$id" "filtered.ordered.ref.fasta"; then
                sed -i "s/\$id/\$label/g" "filtered.ordered.ref.fasta"
            fi
        done < "$ref_seq_list"
        """
}

process GET_FASTA_LEN {
    tag "${target}.on.${reference}"
    container "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"
    
    input:        
        tuple val(target), val(reference), path(filtered_ordered_target_fasta), path(filtered_ordered_ref_fasta)
    
    output:
        tuple val("${target}.on.${reference}"), path("target.seq.lengths"), path("ref.seq.lengths"), emit: tags_len_files
    
    script:
        """
        awk '/^>/ {if (seqlen){print seqlen}; printf \$1" " ;seqlen=0;next; } { seqlen += length(\$0)}END{print seqlen}' $filtered_ordered_target_fasta \
        | sed 's/>//1' \
        | awk '{print \$1, \$2}' OFS='\\t' \
        > target.seq.lengths

        awk '/^>/ {if (seqlen){print seqlen}; printf \$1" " ;seqlen=0;next; } { seqlen += length(\$0)}END{print seqlen}' $filtered_ordered_ref_fasta \
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

    input:
        tuple val(target_on_ref), path(coords_file), path(report_file)
    
    output:
        tuple val(target_on_ref), path("*.1coords.bundle.txt")
    
    script:
        """
        cat $coords_file | awk '{print \$12,\$1,\$2,\$13,\$3,\$4}' OFS="\\t" > "${target_on_ref}.1coords.links.txt"

        /usr/share/circos/tools/bundlelinks/bin/bundlelinks \
        -links "${target_on_ref}.1coords.links.txt" \
        1>"${target_on_ref}.1coords.bundle.txt" \
        2>bundlelinks.err

        add_color_2_circos_bundle_file.pl \
        -i="${target_on_ref}.1coords.bundle.txt" \
        -o="${target_on_ref}.1coords.bundle.coloured.txt"
        """
}

process ADD_COLOUR_TO_BUNDLE_LINKS {
    tag "${target_on_ref}"
    publishDir "${params.outdir.main}/synteny/bundlelinks", mode: 'copy'

    input:
        tuple val(target_on_ref), path(bundle_links)
    
    output:
        tuple val(target_on_ref), path("*.1coords.bundle.coloured.txt")
    
    script:
        """
        add_color_2_circos_bundle_file.pl \
        -i="${bundle_links}" \
        -o="${target_on_ref}.1coords.bundle.coloured.txt"
        """
}

process SPLIT_BUNDLE_FILE_BY_TARGET_SEQS {
    tag "${target_on_ref}"

    input:
        tuple val(target_on_ref), path(coloured_bundle_links)
    
    output:
        tuple val(target_on_ref), path("*.split.bundle.txt")
    
    script:
        """
        target_seqs=(\$(awk '{print \$4}' $coloured_bundle_links | sort | uniq))

        for i in "\${!target_seqs[@]}"
        do
            target_seq=\${target_seqs[\$i]}
            awk -v seq="\$target_seq" '\$4==seq {print \$0}' $coloured_bundle_links > "${target_on_ref}.\${target_seq}.split.bundle.txt"
        done
        """
}

process GENERATE_KARYOTYPE {
    tag "${target_on_ref}.${seq_tag}"

    input:
        tuple val(target_on_ref), val(seq_tag), path(split_bundle_file), path(target_seq_len), path(ref_seq_len)
    
    output:
        tuple val("${target_on_ref}.${seq_tag}"), path("*.karyotype")
    
    script:
        """
        ref_seqs=(\$(awk '{print \$1}' $split_bundle_file | sort | uniq))
        tmp_file=\$(mktemp)
        printf '%s\\n' "\${ref_seqs[@]}" > "\$tmp_file"
        grep "$seq_tag" $target_seq_len > filtered.target.seq.len
        grep -f "\$tmp_file" $ref_seq_len > filtered.ref.seq.len

        paste -d "\\n" filtered.target.seq.len filtered.ref.seq.len > merged.seq.lengths
        sed -i '/^\$/d' merged.seq.lengths
        cat merged.seq.lengths | awk '{print "chr -",\$1,\$1,"0",\$2-1,(\$1=="$seq_tag"?"red":"blue")}' OFS="\t" > "${target_on_ref}.${seq_tag}.karyotype"

        rm "\$tmp_file"
        """
}

process CIRCOS {
    tag "${target_on_ref_seq}"
    container "docker://gallvp/circos-tools:0.23-1" 
    publishDir "${params.outdir.main}/synteny", mode: 'copy'

    input:
        tuple val(target_on_ref_seq), val(karyotype), path(bundle_file)
    
    output:
        path "*.svg", emit: svg_file
        path "*.png", emit: png_file
    
    script:
        """
        cat $bundle_file | awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' OFS="\\t" > bundled.links.tsv

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

        mv circos.svg "${target_on_ref_seq}.svg"
        mv circos.png "${target_on_ref_seq}.png"
        """
}