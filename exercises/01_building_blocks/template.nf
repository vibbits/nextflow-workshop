#!/usr/bin/env nextflow

// Activate the DSL2 - always include at beginning of Nextflow script.
nextflow.enable.dsl=2

// Create the channels
strings_ch = Channel.from('This', 'is', 'a', 'channel')

// Inspect a channels contents with the operator .view()
strings_ch.view()