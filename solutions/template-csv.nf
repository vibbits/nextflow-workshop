#!/usr/bin/env nextflow

// Activate the DSL2 - always include at beginning of Nextflow script.
nextflow.enable.dsl=2

// Create the channels
samples_ch = Channel
                .fromPath('input.csv')
                .splitCsv(header:true)

// Inspect a channels contents with the operator .view()
samples_ch.view{ row -> tuple(row.sampleId, file(row.forward_read), file(row.reverse_read)) }