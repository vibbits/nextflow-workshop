#!/usr/bin/env nextflow

process CUTADAPT {
    // DIRECTIVES: set the docker container, the directory to output to, and a tag to follow along which sample is currently being processed
    container 'kfdrc/cutadapt:latest'
    publishDir 'outputs/trimmed'
    tag "$sample"

    input:
    tuple val( sample ), path( pe_reads )
    val(fw_primer)
    val(rv_primer)
    val(fw_rc_primer)
    val(rv_rc_primer)

    output:
    tuple val( sample ), file( '*_trimmed.fastq' )

    script:
    """
    cutadapt -a ^${fw_primer}...${rv_rc_primer} -A ^${rv_primer}...${fw_rc_primer} --quality-cutoff 28 --minimum-length 30 --max-n 0 --output ${sample}_R1_trimmed.fastq --paired-output ${sample}_R2_trimmed.fastq ${pe_reads[0]} ${pe_reads[1]}
    """
}