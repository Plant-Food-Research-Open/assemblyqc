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

    validateLAIMonoploidSeqs(params)
}

def validateLAIMonoploidSeqs(params) {

    def listOfMonoploidSeqs = params["lai"]["monoploid_seqs"]
    def listOfFastaTuples   = params["target_assemblies"]

    if (listOfMonoploidSeqs.isEmpty()) {
        return
    }

    if (isNotListOfLists(listOfMonoploidSeqs, 2)) {
        error 'Error: lai::monoploid_seqs must be a list of sublists, with each sublist containing 2 elements'
    }

    def fastaTags = listOfFastaTuples.collect { it[0] }
    def monoSeqTags = listOfFastaTuples.collect { it[0] }

    monoSeqTags.each {
        if(!fastaTags.contains(it)) {
            error "Error: $it in lai::monoploid_seqs does not have a corresponding tag in target_assemblies"
        }
    }

    listOfMonoploidSeqs.each {
        validateMonoSeqs(it[1])
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

    xrefTags.each {
        if (fastaTags.contains(it)) {
            error "Error: Tag $it in synteny::xref_assemblies is already included in target_assemblies"
        }
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

def validateMonoSeqs(monoSeqsPath) {
    def monoSeqsFile = file(monoSeqsPath, checkIfExists: true)

    def lines = monoSeqsFile.readLines()
    if (lines.isEmpty()) {
        error "${monoSeqsPath} is empty. It should be a single column text file"
    }

    lines.each { line ->
        def literals = line.split()
        if (literals.size() != 1) {
            error "${monoSeqsPath} should be a single column text file"
        }
    }
}
