import groovy.json.JsonSlurper

def checkHiCParam(paramValue, schema) {
    def jsonSlurper = new JsonSlurper()
    def jsonContent = jsonSlurper.parse ( file ( schema, checkIfExists: true ) )
    def pattern = jsonContent['$defs'].hic_options.properties.hic.pattern
    def match = paramValue ==~ pattern

    return match
}
