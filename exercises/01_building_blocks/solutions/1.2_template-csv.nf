#!/usr/bin/env nextflow

// Create the channels
samples_ch = Channel
                .fromPath('exercises/01_building_blocks/input.csv')
                .splitCsv(header:true)

// Inspect a channels contents with the operator .view()
samples_ch.view()
