#!/usr/bin/env nextflow

fifo = Channel.from(1,2,3,4,5,6,7,8,9,10)

process whosfirst {

    input:
    val x from fifo

    // Using the native execution instead of the more common 'script':
    exec:
    println "This is job number $x"
}

