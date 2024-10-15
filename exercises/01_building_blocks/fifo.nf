#!/usr/bin/env nextflow

// Define the process
process whosfirst {

    input:
    val x 

    // Using the native execution instead of the more common 'script':
    exec:
    println "This is job number $x"
}

workflow {
    // Create channel
    def fifo = Channel.of(1,2,3,4,5,6,7,8,9,10)

    // call the process as a function with channel as its input
    whosfirst(fifo)
}
