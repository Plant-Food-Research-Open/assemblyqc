nextflow.enable.dsl=2

process AGP2_ASSEMBLY {

    container "docker://gallvp/juicebox_scripts:a7ae991"

    input:
        path agp_file

    output:
        path '*.agp.assembly'

    script:
        """
        file_name="$agp_file"
        agp2assembly.py $agp_file "\${file_name%%.*}.agp.assembly"
        """
}