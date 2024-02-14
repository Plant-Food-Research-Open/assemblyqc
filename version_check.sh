#!/usr/bin/env bash

report_version=$(sed -n '/"SELF": /s/"SELF": "v//p' ./bin/assembly_qc_report_943e0fb.py | sed 's/",//' | tr -d '[:space:]')
cff_version=$(sed -n '/^version: /s/version: //p' CITATION.cff | tr -d '[:space:]')

if [[ $report_version != $cff_version ]]; then
    echo 'report_version != cff_version'
    exit 1
fi

# Check README version

readme_citation_text=$(grep -A10 '## Citations' README.md)

if [[ ! $readme_citation_text =~ "($report_version)"  ]]; then
    echo 'Failed to match README version'
    exit 1
fi

# Check CHANGELOG version

grep "Version $report_version" CHANGELOG.md >/dev/null \
    || (echo 'Failed to match CHANGELOG version'; exit 1)
