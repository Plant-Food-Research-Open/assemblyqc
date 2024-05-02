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

workflow FASTA_SYNTENY {
    take:
    ch_fasta                            // Channel: [ tag, fa ]
    ch_labels                           // Channel: [ tag, txt ]
    ch_xref_fasta_labels                // Channel: [ tag2, fa, txt ]
    between_input_assemblies            // val(true|false)
    many_to_many_align                  // val(true|false)
    max_gap                             // val(Integer)
    min_bundle_size                     // val(Integer)
    plot_1_vs_all                       // val(true|false)
    color_by_contig                     // val(true|false)
    plot_type                           // val(linear|circos|both)

    main:
    ch_versions                         = Channel.empty()

    ch_fasta_labels                     = ch_fasta
                                        | join(
                                            ch_labels
                                        )

    ch_input_combinations               = ! between_input_assemblies
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

    ch_xref_ungz_fa_labels              = GUNZIP_FASTA.out.gunzip
                                        | mix(
                                            ch_xref_fa_branch.rest
                                        )
                                        | map { meta, fa -> [ meta.id, fa ] }
                                        | join(
                                            ch_xref_fasta_labels
                                        )
                                        | map { tag, fa, input_fa, seq_list ->
                                            [ tag, fa, seq_list ]
                                        }

    ch_all_combinations                 = ch_input_combinations
                                        | mix(
                                            ch_fasta_labels
                                            | combine(
                                                ch_xref_ungz_fa_labels
                                            )
                                        )

    ch_all_combination_labels           = ch_all_combinations
                                        | map { target_tag, target_fa, target_txt, xref_tag, xref_fa, xref_txt ->
                                            [ "${target_tag}.on.${xref_tag}", target_txt, xref_txt ]
                                        }

    ch_versions                         = ch_versions.mix(GUNZIP_FASTA.out.versions.first())

    // MODULE: FILTERSORTFASTA
    FILTERSORTFASTA ( ch_all_combinations )

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
        many_to_many_align
    )

    ch_versions                         = ch_versions.mix(DNADIFF.out.versions.first())

    // MODULE: BUNDLELINKS
    BUNDLELINKS(
        DNADIFF.out.coords,
        max_gap,
        min_bundle_size
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
                                        | join(ch_all_combination_labels)

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
                                        | join(ch_all_combination_labels)

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
    ch_circos_inputs                    = ( plot_type in [ 'circos', 'both' ] )
                                        ? ch_split_links
                                        | map { target_on_xref, seq_tag, txt ->
                                            [ "${target_on_xref}.${seq_tag}", txt ]
                                        }
                                        | join(GENERATEKARYOTYPE.out.karyotype)
                                        : Channel.empty()
    CIRCOS ( ch_circos_inputs )

    ch_versions                         = ch_versions.mix(CIRCOS.out.versions.first())

    // MODULE: LINEARSYNTENY
    ch_linear_synteny_inputs            = ( plot_type in [ 'linear', 'both' ] )
                                        ? ch_split_links
                                        | map { target_on_xref, seq_tag, txt ->
                                            [ "${target_on_xref}.${seq_tag}", txt ]
                                        }
                                        | join(GENERATEKARYOTYPE.out.karyotype_ref)
                                        | join(GENERATEKARYOTYPE.out.karyotype_target)
                                        : Channel.empty()

    LINEARSYNTENY ( ch_linear_synteny_inputs )

    ch_versions                         = ch_versions.mix(LINEARSYNTENY.out.versions.first())

    emit:
    png                                 = CIRCOS.out.png_file
    html                                = LINEARSYNTENY.out.html
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
