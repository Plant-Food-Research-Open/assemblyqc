nextflow.enable.dsl=2

import groovy.json.JsonOutput

def jsonifyParams(params) {
    return JsonOutput.toJson(params)
}

def validateParams(params) {
    validateFastaTags(params)
    validateGff3Tags(params)
    validateGff3FastaCorrespondence(params)

    validateBUSCOParameters(params)
    validateLAIParameters(params)
    validateSyntenyParameters(params)
}

def validateFastaTags(params) {
    def listOfFastaTuples   = params["target_assemblies"]

    if (isNotListOfLists(listOfFastaTuples, 2)) {
        error 'Error: target_assemblies must be a list of sublists, with each sublist containing 2 elements'
    }

    def fastaTags = listOfFastaTuples.collect { it[0] }

    fastaTags.each {
        if (!(it =~ /^\w+$/)) {
            error "Error: $it is not a valid tag in target_assemblies"
        }
    }

    if (fastaTags.size() != (fastaTags as Set).size()) {
        error "All the tags in target_assemblies should be unique"
    }
}

def validateGff3Tags(params) {
    def listOfGff3Tuples   = params["assembly_gff3"]

    if (listOfGff3Tuples == null) {
        return
    }

    if (listOfGff3Tuples.isEmpty()) {
        return
    }

    if (isNotListOfLists(listOfGff3Tuples, 2)) {
        error 'Error: assembly_gff3 must be a list of sublists, with each sublist containing 2 elements'
    }

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

    def fastaTags = listOfFastaTuples.collect { it[0] }
    def gff3Tags = listOfGff3Tuples.collect { it[0] }

    gff3Tags.each {
        if(!fastaTags.contains(it)) {
            error "Error: $it in assembly_gff3 does not have a corresponding tag in target_assemblies"
        }
    }
}

def validateBUSCOParameters(params) {

    if (params["busco"]["skip"] == 1) {
        return
    }

    listOfBUSCOLineages = params["busco"]["lineage_datasets"]

    if (!(listOfBUSCOLineages instanceof List)) {
        error 'Error: busco::lineage_datasets must be a list of lineage(s)'
    }
}

def validateLAIParameters(params) {

    if (params["lai"]["skip"] == 1) {
        return
    }

    def listOfPassLists     = params["lai"]["pass_list"]
    def listOfOutFiles      = params["lai"]["out_file"]
    def listOfFastaTuples   = params["target_assemblies"]

    if (listOfPassLists == null && listOfOutFiles == null) {
        return
    }

    if (isNotListOfLists(listOfPassLists, 2)) {
        error 'Error: lai::pass_list must be null or a list of sublists, with each sublist containing 2 elements'
    }

    if (isNotListOfLists(listOfOutFiles, 2)) {
        error 'Error: lai::out_file must be null or a list of sublists, with each sublist containing 2 elements'
    }

    if (listOfPassLists.size() != listOfOutFiles.size()) {
        error "Error: The number of elements in lai::pass_list and lai::out_file should be equal"
    }

    if (listOfPassLists.size() != listOfFastaTuples.size()) {
        error "Error: The number of elements in lai::pass_list, lai::out_file and target_assemblies should be equal"
    }

    def passListTags        = listOfPassLists.collect { it[0] }
    def outFileTags         = listOfOutFiles.collect { it[0] }
    def fastaTags           = listOfFastaTuples.collect { it[0] }

    if (!(passListTags.containsAll(fastaTags) && outFileTags.containsAll(fastaTags))) {
        error "Error: The tags in lai::pass_list and lai::out_file should match the tags in target_assemblies"
    }
}

def validateSyntenyParameters(params) {
    if (params["synteny"]["skip"] == 1) {
        return
    }

    def listOfFastaTuples       = params["target_assemblies"]
    def listOfTargetSeqLists    = params["synteny"]["assembly_seq_list"]

    if (isNotListOfLists(listOfTargetSeqLists, 2)) {
        error 'Error: synteny::assembly_seq_list must be a list of sublists, with each sublist containing 2 elements'
    }

    if (listOfTargetSeqLists.size() != listOfFastaTuples.size()) {
        error "Error: The number of elements in synteny::assembly_seq_list and target_assemblies should be equal"
    }
    
    def fastaTags           = listOfFastaTuples.collect { it[0] }
    def seqListTags         = listOfTargetSeqLists.collect { it[0] }

    if (!seqListTags.containsAll(fastaTags)) {
        error "Error: The tags in synteny::assembly_seq_list should match the tags in target_assemblies"
    }

    listOfTargetSeqLists.each {
        validateSeqList(it[1])
    }

    def listOfXRefAssemblies    = params["synteny"]["xref_assemblies"]

    if (listOfXRefAssemblies == null) {
        return
    }

    if (listOfXRefAssemblies.isEmpty()) {
        return
    }

    if (isNotListOfLists(listOfXRefAssemblies, 3)) {
        error 'Error: synteny::xref_assemblies must be a list of sublists, with each sublist containing 3 elements'
    }

    def xrefTags = listOfXRefAssemblies.collect { it[0] }
    
    xrefTags.each {
        if (!(it =~ /^\w+$/)) {
            error "Error: $it is not a valid tag in synteny::xref_assemblies"
        }
    }

    if (xrefTags.size() != (xrefTags as Set).size()) {
        error "All the tags in synteny::xref_assemblies should be unique"
    }

    listOfXRefAssemblies.each {
        validateSeqList(it[2])
    }
}

def isNotListOfLists(thisOne, subListSize) {
    return (!(thisOne instanceof List) || thisOne.isEmpty() || thisOne.any { !(it instanceof List) || it.size() != subListSize })
}

def validateSeqList(seqListPath) {
    def seqListFile = file(seqListPath, checkIfExists: true)

    def lines = seqListFile.readLines()
    if (lines.isEmpty()) {
        error "${seqListPath} is empty. It should be a 2 column tab-delimited file"
    }

    lines.each { line ->
        def columns = line.split("\t")
        if (columns.size() != 2) {
            error "${seqListPath} should be a 2 column tab-delimited file"
        }
    }

    def column1Set = new HashSet<>()
    def column2Set = new HashSet<>()
    def hasUniqueElements = lines.every { line ->
        def columns = line.split("\t")
        def element1 = columns[0]
        def element2 = columns[1]
        column1Set.add(element1) && column2Set.add(element2)
    }

    if (!hasUniqueElements) {
        error "${seqListPath} contains duplicate elements in one or both columns"
    }
}