#!/usr/bin/env bash

module_authors=$(find ./modules -name meta.yml | xargs -I {} grep -A20 'authors:' {} | grep '\- ' | tr -d '"' | tr '[:upper:]' '[:lower:]' | awk '{print $2}')
workflow_authors=$(find ./subworkflows -name meta.yml | xargs -I {} grep -A20 'authors:' {} | grep '\- ' | tr -d '"' | tr '[:upper:]' '[:lower:]' | awk '{print $2}')
echo -e "${module_authors}\n${workflow_authors}" | sort -V | uniq | sed -n 's|@\(.*\)|<a href="https://github.com/\1"><img src="https://github.com/\1.png" width="50" height="50"></a>|p'
