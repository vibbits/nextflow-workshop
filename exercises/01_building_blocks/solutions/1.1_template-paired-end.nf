#!/usr/bin/env nextflow

// Create the channels
paired_ch = Channel.fromFilePairs('./data/*{1,2}.fq.gz')

// Inspect a channels contents with the operator .view()
paired_ch.view()