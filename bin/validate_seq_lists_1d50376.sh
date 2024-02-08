#!/usr/bin/env bash

seqFileA=$1
seqFileB=$2

linesFileA=()
linesFileB=()

while IFS= read -r line; do
    linesFileA+=("$line")
    columns=($line)
    if [[ ${#columns[@]} -ne 2 ]]; then
        echo "Error: Sequence file $(basename "$seqFileA") does not have exactly two columns." >&2
        exit 1
    fi
done < "$seqFileA"

while IFS= read -r line; do
    linesFileB+=("$line")
    columns=($line)
    if [[ ${#columns[@]} -ne 2 ]]; then
        echo "Error: Sequence file $(basename "$seqFileB") does not have exactly two columns." >&2
        exit 1
    fi
done < "$seqFileB"

outputLines=("${linesFileA[@]}" "${linesFileB[@]}")

secondColumn=()
for line in "${outputLines[@]}"; do
    columns=($line)
    secondColumn+=("${columns[1]}")
done

uniqueSecondColumn=($(echo "${secondColumn[@]}" | tr ' ' '\n' | sort -u))
if [[ ${#secondColumn[@]} -ne ${#uniqueSecondColumn[@]} ]]; then
    echo "Error: Duplicate sequence labels detected in second column for pair: $(basename "$seqFileA"), $(basename "$seqFileB")" >&2
    exit 1
fi
