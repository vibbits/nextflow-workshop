#!/usr/bin/env nextflow

// this is an example of how to create and update a meta-map
params.name = "Kobe"
params.age = 26


meta = [:]
meta.name = params.name
meta.age = params.age
meta.number = 44
println meta


process testprocess {
    input:
    val(meta)

    output:
    stdout

    """
    echo ${meta.name}
    """

}


workflow {
    testprocess(meta)
    testprocess.out.view()
}