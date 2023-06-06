#!/usr/bin/env bash

## Inputs
fasta_file="$1"
gff3_file="$2"

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
    echo "Failed to validate gff3 file for: $tag_label"
    echo "Fasta file: $fasta_file"
    echo "Gff3 file: $gff3_file"
    echo "GFF3 file contains identifiers not present in FASTA:"
    echo "$missing_identifiers"
    exit 1
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
        echo "Failed to validate gff3 file for: $tag_label"
        echo "Fasta file: $fasta_file"
        echo "Gff3 file: $gff3_file"
        echo "Coordinates exceed sequence length in GFF3 file:"
        echo "Sequence: $seqname"
        echo "Sequence length: $seq_length"
        echo "Start: $start"
        echo "End: $end"
        exit 1
    fi
done < "$gff3_file" 3< "sequence_lengths.txt"