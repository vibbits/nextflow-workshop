#!/usr/bin/env nextflow

params.reads = "${launchDir}/data/*{1,2}.fq.gz"

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromFilePairs( params.reads, checkIfExists:true )
    .view()
    
process fastqc {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'

    input:
    tuple val(sample), path(read)  
    
    script:
    """
    fastqc ${read}
    """
}

workflow {
    fastqc(reads_ch)
}