#!/usr/bin/env nextflow

process DADA2 {
    // DIRECTIVES: set the docker container, the directory to output to, and the amount of cpus to allocate to this process
    container 'golob/dada2:1.14.1.ub.1804'
    publishDir "${params.outdir}/dada2", mode: 'copy', overwrite: true
    cpus 4

    input:
    path(reads)

    output:
    path('counts_matrix.csv')
    path('dendrogram.png')

    script:
    """
    reads2counts.r ${reads}
    """
}


process DADA2_ALTERNATIVE {
    // DIRECTIVES: set the docker container, the directory to output to, and the amount of cpus to allocate to this process
    container 'golob/dada2:1.14.1.ub.1804'
    publishDir "${params.outdir}/dada2_alternative", mode: 'copy', overwrite: true
    cpus 4

    input:
    path(reads)
    path(script)

    output:
    path('counts_matrix.csv')
    path('dendrogram.png')

    script:
    """
    Rscript ${script} ${reads}
    """
}