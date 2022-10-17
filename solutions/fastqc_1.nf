#!/usr/bin/env nextflow

params.reads = "$launchDir/../../data/*.fq.gz"

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromPath( params.reads )
    .view()
    
process fastqc {

    input:
    file read  
    
    script:
    """
    fastqc ${read}
    """
}

workflow{
    fastqc(reads_ch)
}