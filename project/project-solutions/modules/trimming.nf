#!/usr/bin/env nextflow

process CUTADAPT {
    // DIRECTIVES: set the docker container, the directory to output to, and a tag to follow along which sample is currently being processed
    container 'kfdrc/cutadapt:latest'
    publishDir "${params.outdir}/cutadapt", mode: 'copy', overwrite: true
    tag "${sample}"

    input:
    tuple val(sample), path(pe_reads)

    output:
    tuple val(sample), file('*_trimmed.fastq')

    script:
    """
    cutadapt -a ^${params.fw_primer}...${params.rv_primer_rev_comp} \
        -A ^${params.rv_primer}...${params.fw_primer_rev_comp} \
        --quality-cutoff 28 \
        --minimum-length 30 \
        --max-n 0 \
        --output ${sample}_R1_trimmed.fastq \
        --paired-output ${sample}_R2_trimmed.fastq \
        ${pe_reads[0]} ${pe_reads[1]}
    """
}