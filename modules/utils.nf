nextflow.enable.dsl=2

import groovy.json.JsonOutput

def jsonifyParams(params) {
    return JsonOutput.toJson(params)
}

def validateParams(params) {
    validateFastaTags(params)
    validateGff3Tags(params)
    validateGff3FastaCorrespondence(params)
}

def validateFastaTags(params) {
    def listOfFastaTuples   = params["target_assemblies"]

    def fastaTags = listOfFastaTuples.collect { it[0] }

    fastaTags.each {
        if (!(it =~ /^\w+$/)) {
            error "Error: $it is not a valid tag in target_assemblies"
        }
    }
}

def validateGff3Tags(params) {
    def listOfGff3Tuples   = params["assembly_gff3"]

    def gff3Tags = listOfGff3Tuples.collect { it[0] }

    gff3Tags.each {
        if (!(it =~ /^\w+$/)) {
            error "Error: $it is not a valid tag in assembly_gff3"
        }
    }
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

    gff3Tags.each {
        if(!fastaTags.contains(it)) {
            error "Error: $it in assembly_gff3 does not have a corresponding tag in target_assemblies"
        }
    }
}