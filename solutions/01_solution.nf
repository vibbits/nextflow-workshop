#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_csv = 'exercises/input.csv'

samples_ch = Channel
                .fromPath(params.input_csv)
                .splitCsv(header:true)
                .map{ row-> tuple(row.sampleId, file(row.forward_read), file(row.reverse_read)) }

process split_csv {
    input:
    tuple val(sampleId), file(read1), file(read2)  

    script:
    """
    echo your_command --sample $sampleId --reads $read1 $read2 
    """
}

workflow{
    samples_ch.view()
    split_csv(samples_ch)
}