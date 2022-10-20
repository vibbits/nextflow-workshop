#!/usr/bin/env nextflow

params.reads = "${launchDir}/data/*{1,2}.fq.gz"
params.outdir = "${launchDir}/results"

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromFilePairs( params.reads, checkIfExists:true )
    .view()
    
process fastqc {
    publishDir "${params.outdir}/quality-control-${sample}/", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/fastqc:0.11.9--0'


    input:
    tuple val(sample), path(read)  
    
    output:
    path("*_fastqc.{zip,html}") 
    
    script:
    """
    fastqc ${read}
    """
}

workflow {
    fastqc(reads_ch)
}