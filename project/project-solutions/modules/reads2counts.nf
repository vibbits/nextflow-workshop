#!/usr/bin/env nextflow

process DADA2 {
    // DIRECTIVES: set the docker container, the directory to output to, and the amount of cpus to allocate to this process
    container 'golob/dada2:1.14.1.ub.1804'
    publishDir 'outputs'
    cpus 4

    input:
    path( reads )
    path( script )

    output:
    path( 'counts_matrix.csv' )
    path( 'dendrogram.png' )

    script:
    def reads_arg = reads.join(' ')
    """
    Rscript $script $reads_arg
    """
}