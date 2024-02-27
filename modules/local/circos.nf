process CIRCOS {
    tag "${target_on_ref_seq}"
    label 'process_single'

    container "docker.io/gallvp/circos-tools:v0.23-1_ps"

    input:
    tuple val(target_on_ref_seq), path(karyotype), path(bundle_file)

    output:
    path "*.svg", emit: svg_file
    path "*.png", emit: png_file
    path "bundled.links.tsv", emit: bundled_links_tsv
    path "circos.conf", emit: circos_conf
    path "karyotype.tsv", emit: karyotype_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

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

    cat <<-END_CONF > circos.conf
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
    END_CONF

    circos

    mv circos.svg "${target_on_ref_seq}.svg"
    mv circos.png "${target_on_ref_seq}.png"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circos: \$(circos -version | awk '{print \$2}' FS='|' | tr -d '[:space:]')
        perl: \$(circos -version | awk '{print \$4}' FS='|' | tr -d '[:space:]Perl')
    END_VERSIONS
    """

    stub:
    """
    touch ${target_on_ref_seq}.svg
    touch ${target_on_ref_seq}.png

    touch bundled.links.tsv
    touch circos.conf
    touch karyotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circos: \$(circos -version | awk '{print \$2}' FS='|' | tr -d '[:space:]')
        perl: \$(circos -version | awk '{print \$4}' FS='|' | tr -d '[:space:]Perl')
    END_VERSIONS
    """
}
