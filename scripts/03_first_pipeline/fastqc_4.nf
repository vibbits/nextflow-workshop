#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "$launchDir/../../data/*{1,2}.fq.gz"
params.outdir = "$launchDir/results"

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromFilePairs( params.reads, checkIfExists:true )
    .view()
    
process fastqc {
    publishDir "$params.outdir/quality-control-$sample/", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read)  
    
    output:
    path("*_fastqc.{zip,html}") 
    
    script:
    """
    fastqc ${read}
    """
}
//mkdir -p $params.outdir/quality-control-$sample/

workflow{
    fastqc(reads_ch)
}