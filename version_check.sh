#!/usr/bin/env bash

config_version=$(sed -n "/^\s*version\s*=\s*'/s/version//p" nextflow.config | tr -d "=[:space:]'")
cff_version=$(sed -n '/^version: /s/version: //p' CITATION.cff | tr -d '[:space:]')

if [[ $config_version != $cff_version ]]; then
    echo 'config_version != cff_version'
    exit 1
fi

# Check CHANGELOG version

head -10 CHANGELOG.md | grep "## v$config_version - " >/dev/null \
    || (echo 'Failed to match CHANGELOG version'; exit 1)

# Check README version

grep "assembly quality ($config_version)" README.md >/dev/null \
    || (echo 'Failed to match README version'; exit 1)
