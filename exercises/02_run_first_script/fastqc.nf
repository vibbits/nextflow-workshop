#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "$launchDir/../../data/*.fq.gz" // $projectDir is another interesting implicit variable


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