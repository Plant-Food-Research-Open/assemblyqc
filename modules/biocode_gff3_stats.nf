nextflow.enable.dsl=2

process BIOCODE_GFF3_STATS {
    tag "${tag_label}"
    label "process_single"
    
    container "quay.io/biocontainers/biocode:0.10.0--pyhdfd78af_0"
    publishDir "${params.outdir.main}/biocode_gff3_stats", mode: 'copy'

    input:
        tuple val(tag_label), path(gff3_file)

    output:
        path "${tag_label}_stats.csv"

    script:
        """
        valid_regions=("gene" "mrna" "cds" "exon")
        validity=true

        while IFS=\$'\\t' read -r _ _ region _
        do
        region=\$(echo "\$region" | tr '[:upper:]' '[:lower:]')
        if [[ ! " \${valid_regions[@]} " =~ " \$region " ]]; then
            validity=false
            break
        fi
        done < "$gff3_file"

        if [ "\$validity" = true ]; then
            report_gff3_statistics.py --input_file "$gff3_file" > "${tag_label}_stats.csv"
        else
            echo "Failed to compute statistics. This module expects a 3-level gff3 file with following levels: gene/mRNA/exon,CDS" \
            > "${tag_label}_stats.csv"
        fi
        """ 
}