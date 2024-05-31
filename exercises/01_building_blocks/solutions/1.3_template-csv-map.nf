#!/usr/bin/env nextflow

// Create the channels
samples_ch = Channel
                .fromPath('exercises/01_building_blocks/input.csv')
                .splitCsv(header:true)

// Inspect a channels contents with the operator .view()
samples_ch.view()
samples_ch.view{ row -> tuple(row.sampleId, file(row.forward_read), file(row.reverse_read)) }

// read1_ch = samples_ch.map{ row -> tuple(row.sampleId, row.forward_read) }

// By using the file() function, you create a file object from the path strings in the csv file:
read1_ch = samples_ch.map{ row -> tuple(row.sampleId, file(row.forward_read)) }

read2_ch = samples_ch.map{ row -> tuple(row.sampleId, file(row.reverse_read)) }

read1_ch.view()
read2_ch.view()
