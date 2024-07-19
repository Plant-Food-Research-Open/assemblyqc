include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/gunzip/main'
include { FILTERSORTFASTA               } from '../../modules/local/filtersortfasta'
include { MUMMER                        } from '../../modules/local/mummer'
include { GETFASTALENGTH                } from '../../modules/local/getfastalength'
include { DNADIFF                       } from '../../modules/local/dnadiff'
include { BUNDLELINKS                   } from '../../modules/local/bundlelinks'
include { COLOURBUNDLELINKS             } from '../../modules/local/colourbundlelinks'
include { RELABELBUNDLELINKS            } from '../../modules/local/relabelbundlelinks'
include { SPLITBUNDLEFILE               } from '../../modules/local/splitbundlefile'
include { RELABELFASTALENGTH            } from '../../modules/local/relabelfastalength'
include { GENERATEKARYOTYPE             } from '../../modules/local/generatekaryotype'
include { CIRCOS                        } from '../../modules/local/circos'
include { LINEARSYNTENY                 } from '../../modules/local/linearsynteny'
include { CUSTOM_RELABELFASTA           } from '../../modules/pfr/custom/relabelfasta/main'
include { MINIMAP2_ALIGN                } from '../../modules/nf-core/minimap2/align/main'
include { SYRI                          } from '../../modules/pfr/syri/main'
include { PLOTSR                        } from '../../modules/pfr/plotsr/main'

workflow FASTA_SYNTENY {
    take:
    ch_fasta                            // Channel: [ tag, fa ]
    ch_labels                           // Channel: [ tag, txt ]
    ch_xref_fasta_labels                // Channel: [ tag2, fa, txt ]
    between_input_assemblies            // val(true|false)
    mummer_m2m_align                    // val(true|false)
    mummer_max_gap                      // val(Integer)
    mummer_min_bundle_size              // val(Integer)
    plot_1_vs_all                       // val(true|false)
    color_by_contig                     // val(true|false)
    mummer_plot_type                    // val(linear|circos|both)
    mummer_skip                         // val(true|false)
    plotsr_seq_label                    // val(String)
    plotsr_skip                         // val(true|false)
    plotsr_assembly_order               // val(String)

    main:
    ch_versions                         = Channel.empty()

    ch_fasta_labels                     = ch_fasta
                                        | join(
                                            ch_labels
                                        )

    ch_input_combination                = ! between_input_assemblies
                                        ? Channel.empty()
                                        : ch_fasta_labels
                                        | map { [it] }
                                        | collect
                                        | map { getUniqueWithinCombinations(it) }
                                        | flatten
                                        | buffer(size:6)

    ch_xref_fa_branch                   = ch_xref_fasta_labels
                                        | map { tag, fa, txt ->
                                            [ [ id: tag ], fa ]
                                        }
                                        | branch { meta, fa ->
                                            gz: "$fa".endsWith(".gz")
                                            rest: !"$fa".endsWith(".gz")
                                        }

    // MODULE: GUNZIP_FASTA
    GUNZIP_FASTA ( ch_xref_fa_branch.gz )

    ch_xref_ungz                        = GUNZIP_FASTA.out.gunzip
                                        | mix(
                                            ch_xref_fa_branch.rest
                                        )
                                        | map { meta, fa -> [ meta.id, fa ] }

    ch_xref_ungz_fa_labels              = ch_xref_ungz
                                        | join(
                                            ch_xref_fasta_labels
                                        )
                                        | map { tag, fa, input_fa, seq_list ->
                                            [ tag, fa, seq_list ]
                                        }

    ch_combination                      = mummer_skip
                                        ? Channel.empty()
                                        : ch_input_combination
                                        | mix(
                                            ch_fasta_labels
                                            | combine(
                                                ch_xref_ungz_fa_labels
                                            )
                                        )

    ch_combination_labels               = ch_combination
                                        | map { target_tag, target_fa, target_txt, xref_tag, xref_fa, xref_txt ->
                                            [ "${target_tag}.on.${xref_tag}", target_txt, xref_txt ]
                                        }

    ch_versions                         = ch_versions.mix(GUNZIP_FASTA.out.versions.first())

    // MODULE: FILTERSORTFASTA
    FILTERSORTFASTA ( ch_combination )

    ch_versions                         = ch_versions.mix(FILTERSORTFASTA.out.versions.first())

    // MODULE: MUMMER
    MUMMER ( FILTERSORTFASTA.out.fasta )

    ch_versions                         = ch_versions.mix(MUMMER.out.versions.first())

    // MODULE: GETFASTALENGTH
    GETFASTALENGTH ( FILTERSORTFASTA.out.fasta )

    ch_versions                         = ch_versions.mix(GETFASTALENGTH.out.versions.first())

    // MODULE: DNADIFF
    ch_dnadiff_inputs                   = FILTERSORTFASTA.out.fasta
                                        | map { target, reference, target_fasta, ref_fasta ->
                                            [ "${target}.on.${reference}", target_fasta, ref_fasta ]
                                        }
                                        | join(
                                            MUMMER.out.delta
                                        )
    DNADIFF(
        ch_dnadiff_inputs,
        mummer_m2m_align
    )

    ch_versions                         = ch_versions.mix(DNADIFF.out.versions.first())

    // MODULE: BUNDLELINKS
    BUNDLELINKS(
        DNADIFF.out.coords,
        mummer_max_gap,
        mummer_min_bundle_size
    )

    ch_versions                         = ch_versions.mix(BUNDLELINKS.out.versions.first())

    // MODULE: COLOURBUNDLELINKS
    COLOURBUNDLELINKS(
        BUNDLELINKS.out.links,
        color_by_contig
    )

    ch_coloured_links                   = COLOURBUNDLELINKS.out.coloured_links
    ch_versions                         = ch_versions.mix(COLOURBUNDLELINKS.out.versions.first())

    // MODULE: RELABELBUNDLELINKS
    ch_relabellinks_inputs              = ch_coloured_links
                                        | join(ch_combination_labels)

    RELABELBUNDLELINKS ( ch_relabellinks_inputs )

    ch_versions                         = ch_versions.mix(RELABELBUNDLELINKS.out.versions.first())

    // MODULE: SPLITBUNDLEFILE
    SPLITBUNDLEFILE(
        RELABELBUNDLELINKS.out.relabeled_links,
        plot_1_vs_all
    )

    ch_split_links                      = SPLITBUNDLEFILE.out.split_file
                                        | map { flattenSplitBundles(it) }
                                        | flatten
                                        | buffer(size:3)

    ch_versions                         = ch_versions.mix(SPLITBUNDLEFILE.out.versions.first())

    // MODULE: RELABELFASTALENGTH
    ch_relabelfastalength_inputs        = GETFASTALENGTH.out.length
                                        | join(ch_combination_labels)

    RELABELFASTALENGTH ( ch_relabelfastalength_inputs )

    ch_versions                         = ch_versions.mix(RELABELFASTALENGTH.out.versions.first())

    // MODULE: GENERATEKARYOTYPE
    ch_generate_karyotype_inputs        = RELABELFASTALENGTH.out.relabeled_seq_lengths
                                        | cross(
                                            ch_split_links
                                        )
                                        | map { seq_len_tuple, split_bundle_tuple ->

                                            def target_on_xref      = seq_len_tuple[0]
                                            def seq_tag             = split_bundle_tuple[1]
                                            def split_bundle_file   = split_bundle_tuple[2]
                                            def target_seq_len      = seq_len_tuple[1]
                                            def ref_seq_len         = seq_len_tuple[2]

                                            [ target_on_xref, seq_tag, split_bundle_file, target_seq_len, ref_seq_len ]
                                        }
    GENERATEKARYOTYPE ( ch_generate_karyotype_inputs )

    ch_versions                         = ch_versions.mix(GENERATEKARYOTYPE.out.versions.first())

    // MODULE: CIRCOS
    ch_circos_inputs                    = ( mummer_plot_type in [ 'circos', 'both' ] )
                                        ? ch_split_links
                                        | map { target_on_xref, seq_tag, txt ->
                                            [ "${target_on_xref}.${seq_tag}", txt ]
                                        }
                                        | join(GENERATEKARYOTYPE.out.karyotype)
                                        : Channel.empty()
    CIRCOS ( ch_circos_inputs )

    ch_versions                         = ch_versions.mix(CIRCOS.out.versions.first())

    // MODULE: LINEARSYNTENY
    ch_linear_synteny_inputs            = ( mummer_plot_type in [ 'dotplot', 'both' ] )
                                        ? ch_split_links
                                        | map { target_on_xref, seq_tag, txt ->
                                            [ "${target_on_xref}.${seq_tag}", txt ]
                                        }
                                        | join(GENERATEKARYOTYPE.out.karyotype_ref)
                                        | join(GENERATEKARYOTYPE.out.karyotype_target)
                                        : Channel.empty()

    LINEARSYNTENY ( ch_linear_synteny_inputs )

    ch_versions                         = ch_versions.mix(LINEARSYNTENY.out.versions.first())

    // Create chr label lists
    ch_assembly_labels                  = plotsr_skip
                                        ? Channel.empty()
                                        : ch_fasta_labels
                                        | mix(ch_xref_ungz_fa_labels)

    ch_common_label_count               = ch_assembly_labels
                                        | map { tag, fa, labels ->
                                            labels.readLines().findAll { it != '' }.size()
                                        }
                                        | collect
                                        | map { it.min() }

    ch_plotsr_formatted_labels          = ch_assembly_labels
                                        | combine(ch_common_label_count)
                                        | map { tag, fa, labels, num ->
                                            def label_lines = labels
                                                .readLines()
                                                .collect { it.trim() }
                                                .findAll { it != '' }

                                            label_lines.each { line ->
                                                if ( line.split('\t').size() != 2 ) {
                                                    error "synteny_labels file ${labels.name} for assembly ${tag} is malformed near line ${line}"
                                                }
                                            }

                                            def new_labels = label_lines[0..<num].withIndex().collect { line, index ->
                                                def literals = line.split('\t')
                                                def seq = literals[0]
                                                def label = "${plotsr_seq_label}${index+1}"

                                                "$seq\t$label"
                                            }.join('\n')

                                            [ "${tag}.plotsr.csv", new_labels]
                                        }
                                        | collectFile(storeDir: "${params.outdir}/synteny/plotsr")
                                        | map { labels -> [ labels.baseName.replace('.plotsr', ''), labels ] }

    // MODULE: CUSTOM_RELABELFASTA
    ch_relabel_inputs                   = ch_assembly_labels
                                        | join(ch_plotsr_formatted_labels)
                                        | map { tag, fa, old_l, new_l -> [ tag, fa, new_l ] }
    CUSTOM_RELABELFASTA(
        ch_relabel_inputs.map { tag, fa, labels -> [ [ id: tag ], fa ] },
        ch_relabel_inputs.map { tag, fa, labels -> labels }
    )

    ch_plotsr_assembly                  = CUSTOM_RELABELFASTA.out.fasta
    ch_versions                         = ch_versions.mix(CUSTOM_RELABELFASTA.out.versions.first())

    // MODULE: MINIMAP2_ALIGN
    ch_minimap_inputs                   = ch_plotsr_assembly
                                        | map { [ it ] }
                                        | collect
                                        | map { list ->
                                            if ( plotsr_assembly_order == null ) {
                                                return list.sort(false) { it[0].id.toUpperCase() }
                                            }

                                            def order   = plotsr_assembly_order.tokenize(' ')

                                            if ( order.size() != order.unique().size() ) error "Tags listed by synteny_plotsr_assembly_order should all be unique: $order"

                                            def tags    = list.collect { it[0].id }

                                            def ordered_list = []
                                            order.each { tag ->
                                                if ( ! ( tag in tags ) ) error "Assembly $tag listed in synteny_plotsr_assembly_order could not be found in input or synteny_xref_assemblies: $tags"

                                                ordered_list << ( list.find { it[0].id == tag } )
                                            }
                                            ordered_list
                                        }
                                        | flatMap { list ->
                                            if ( list.size() < 2 ) return null

                                            def new_list = []
                                            list.eachWithIndex { assembly, index -> if ( index > 0 ) { new_list << [ assembly, list[index-1] ] } }
                                            new_list
                                        }
                                        | flatten
                                        | buffer(size: 4)
                                        | map { tmeta, tfa, rmeta, rfa ->
                                            [ [ id: "${tmeta.id}.on.${rmeta.id}" ], tfa, rfa ]
                                        }
    MINIMAP2_ALIGN(
        ch_minimap_inputs.map { meta, tfa, rfa -> [ meta, tfa ] },
        ch_minimap_inputs.map { meta, tfa, rfa -> [ meta, rfa ] },
        true,   // bam_format
        false,  // cigar_paf_format
        false   // cigar_bam
    )

    ch_minimap2_bam                     = MINIMAP2_ALIGN.out.bam
    ch_versions                         = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    // MODULE: SYRI
    ch_syri_inputs                      = ch_minimap2_bam
                                        | join(ch_minimap_inputs)

    SYRI(
        ch_syri_inputs.map { meta, bam, tfa, rfa -> [ meta, bam ] },
        ch_syri_inputs.map { meta, bam, tfa, rfa -> tfa },
        ch_syri_inputs.map { meta, bam, tfa, rfa -> rfa },
        'B' // BAM
    )

    ch_syri                             = SYRI.out.syri
    ch_syri_fail_log                    = SYRI.out.error
    ch_versions                         = ch_versions.mix(SYRI.out.versions.first())

    // MODULE: PLOTSR
    ch_plotsr_inputs                    = ch_syri
                                        | join(ch_syri_inputs)
                                        | map { meta, syri, bam, tfa, rfa ->

                                            [
                                                [ id: 'plotsr' ],
                                                syri,
                                                [ tfa, rfa]
                                            ]
                                        }
                                        | groupTuple
                                        | map { meta, syri, fastas ->
                                            def fasta_list          = fastas.flatten()
                                            def syri_tags           = syri.collect { it.name.replace('syri.out', '').split(/\.on\./).toList() }.flatten().unique()
                                            def proposed_order      = plotsr_assembly_order ? plotsr_assembly_order.tokenize(' ') : syri_tags.sort(false)

                                            def available_tags      = []
                                            proposed_order.each { tag -> if ( tag in syri_tags ) available_tags << tag }

                                            def ordered_fa          = []
                                            available_tags.each { tag -> ordered_fa << ( fasta_list.find { it.baseName == "${tag}.plotsr" } ) }

                                            def ordered_syri_tags   = []
                                            available_tags.eachWithIndex { tag, index -> if ( index > 0 ) { ordered_syri_tags << "${tag}.on.${available_tags[index-1]}" } }

                                            def ordered_syri        = []
                                            ordered_syri_tags.each { tag -> ordered_syri << ( syri.find { it.baseName == "${tag}syri" } ) }

                                            [
                                                meta,
                                                ordered_syri,
                                                ordered_fa,
                                                "#file\tname\n" + ordered_fa.collect { it.baseName.replace('.plotsr', '') }.join('\n')
                                            ]
                                        }

    PLOTSR(
        ch_plotsr_inputs.map { meta, syri, fastas, txt -> [ meta, syri ] },
        ch_plotsr_inputs.map { meta, syri, fastas, txt -> fastas },
        ch_plotsr_inputs.map { meta, syri, fastas, txt -> txt },
        [],
        [],
        [],
        [],
        []
    )

    ch_plotsr_png                       = PLOTSR.out.png
    ch_versions                         = ch_versions.mix(PLOTSR.out.versions.first())

    emit:
    png                                 = CIRCOS.out.png_file
                                        | mix( ch_plotsr_png.map { meta, png -> png } )
    html                                = LINEARSYNTENY.out.html
    syri_fail_log                       = ch_syri_fail_log.map { meta, log -> log }
    plotsr_labels                       = ch_plotsr_formatted_labels.map { tag, labels -> labels }
    versions                            = ch_versions
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
                inputArray[i][1],
                inputArray[i][2],
                inputArray[j][0],
                inputArray[j][1],
                inputArray[j][2]
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
        // This branch should not execut unless the upstream logic is flawed
        error "Error: Failed to parse the sequence tag from file name: ${filePath.getName()}"
    }
}
