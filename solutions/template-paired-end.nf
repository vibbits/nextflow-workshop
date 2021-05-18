#!/usr/bin/env nextflow

// Activate the DSL2 - always include at beginning of Nextflow script.
nextflow.enable.dsl=2

// Create the channels
paired_ch = Channel.fromFilePairs('../data/*{1,2}.fastq')

// Inspect a channels contents with the operator .view()
paired_ch.view()