#!/usr/bin/env bash

# Bump VERSION on edit
VERSION="v1"

gff3_file="!{gff3}"
fasta_file="!{fasta}"
out_prefix="!{prefix}"
task_process="!{task.process}"

# Record versions
cat <<-END_VERSIONS > versions.yml
"${task_process}":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//' )
END_VERSIONS

# Requires
# samtools faidx

## STEP 1
# Check that gff3 has no identifers that are not in fasta (fasta can
# have ids that are not in gff3 since not all assembly units have gff3 records

# Extract identifiers from the GFF3 file
gff3_identifiers=$(grep -v '^#' "$gff3_file" | awk '{print $1}' | sort -u)

# Extract identifiers from the FASTA file
fasta_identifiers=$(grep '^>' "$fasta_file" | awk '{print substr($1, 2)}' | sort -u)

# Compare identifiers and find any that are present in the GFF3 but not in the FASTA
missing_identifiers=$(comm -23 <(echo "$gff3_identifiers") <(echo "$fasta_identifiers"))

# Check if any missing identifiers were found
if [[ -n "$missing_identifiers" ]]; then
    touch "${out_prefix}.error.log"
    echo "Failed to validate gff3 file for: $tag_label"         >> "${out_prefix}.error.log"
    echo "Fasta file: $fasta_file"                              >> "${out_prefix}.error.log"
    echo "Gff3 file: $gff3_file"                                >> "${out_prefix}.error.log"
    echo "GFF3 file contains identifiers not present in FASTA:" >> "${out_prefix}.error.log"
    echo "$missing_identifiers"                                 >> "${out_prefix}.error.log"
    exit 0
fi

## STEP 2
# check that there are no coordiantes in gff3 for any seqid that are
# greater than the seq length of the paretn fasta entry

# Compute sequence lengths using samtools faidx
samtools faidx "$fasta_file" | cut -f 1,2 > sequence_lengths.txt

# Check GFF3 file for coordinates exceeding sequence lengths
while IFS=$'\t' read -r seqname source feature start end score strand frame attributes && \
    read -r seq seq_length <&3; do
    if [[ $start -gt $seq_length || $end -gt $seq_length ]]; then
        touch "${out_prefix}.error.log"
        echo "Failed to validate gff3 file for: $tag_label"     >> "${out_prefix}.error.log"
        echo "Fasta file: $fasta_file"                          >> "${out_prefix}.error.log"
        echo "Gff3 file: $gff3_file"                            >> "${out_prefix}.error.log"
        echo "Coordinates exceed sequence length in GFF3 file:" >> "${out_prefix}.error.log"
        echo "Sequence: $seqname"                               >> "${out_prefix}.error.log"
        echo "Sequence length: $seq_length"                     >> "${out_prefix}.error.log"
        echo "Start: $start"                                    >> "${out_prefix}.error.log"
        echo "End: $end"                                        >> "${out_prefix}.error.log"
        exit 0
    fi
done < "$gff3_file" 3< "sequence_lengths.txt"

touch "${out_prefix}.success.log"
echo "All tests passed..."                                      >>  "${out_prefix}.success.log"
exit 0
