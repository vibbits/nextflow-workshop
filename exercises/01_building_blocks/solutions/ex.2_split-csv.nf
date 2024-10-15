#!/usr/bin/env nextflow

params.input_csv = 'exercises/01_building_blocks/input.csv'

process split_csv {
    input:
    tuple val(sampleId), path(read1), path(read2)  

    script:
    """
    echo your_command --sample $sampleId --reads $read1 $read2 
    """
}

workflow {
    def samples_ch = Channel
                    .fromPath(params.input_csv)
                    .splitCsv(header:true)
                    .map{ row-> tuple(row.sampleId, file(row.forward_read), file(row.reverse_read)) }

    samples_ch.view()
    split_csv(samples_ch)
}