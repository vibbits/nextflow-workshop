#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process FASTQC { 
    // DIRECTIVES: set the docker container, the directory to output to, and a tag to follow along which sample is currently being processed
    container 'biocontainers/fastqc:v0.11.9_cv8'
    publishDir 'outputs/fastqc_logs'
    tag "$sample"

    input:
    tuple val( sample ), path( reads )

    output:
    path( '*_fastqc.{zip,html}' )

    script:
    """
    fastqc -f fastq -q ${reads}
    """
}


process MULTIQC {
    // DIRECTIVES: set the docker container, the directory to output to, and a tag
    container 'ewels/multiqc'
    publishDir 'outputs'
    tag "$state"

    input:
    val( state )
    path( 'fastqc_logs/*' )

    output:
    file( '*_multiqc_report.html' )

    script:
    """
    multiqc --filename ${state}_multiqc_report fastqc_logs
    """
}


// combine the two processes into a subworkflow
workflow check_QC {
    take:
    state
    reads_ch

    main:
    FASTQC( reads_ch )
    MULTIQC( state, FASTQC.out.collect() )
}