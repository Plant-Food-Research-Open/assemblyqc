nextflow.enable.dsl=2

import groovy.json.JsonOutput

def jsonifyParams(params) {
    return JsonOutput.toJson(params)
}

def validateParams(params) {
    validateGff3FastaCorrespondence(params)
}

def validateGff3FastaCorrespondence(params) {
    
    def listOfGff3Tuples    = params["assembly_gff3"]
    def listOfFastaTuples   = params["target_assemblies"]

    if (listOfGff3Tuples == null) {
        return
    }

    if (listOfGff3Tuples.isEmpty()) {
        return
    }

    if (!(listOfGff3Tuples instanceof List) || listOfGff3Tuples.any { !(it instanceof List) || it.size() != 2 }) {
        error 'Error: assembly_gff3 is not configured correctly'
    }

    def fastaTags = listOfFastaTuples.collect { it[0] }
    def gff3Tags = listOfGff3Tuples.collect { it[0] }

    if (!gff3Tags.every { fastaTags.contains(it) }) {
        error 'Error: Not every gff3 tag in assembly_gff3 has a corresponding fasta tag in target_assemblies'
    }

    return
}