#!/usr/bin/env nextflow


// Create the channels
strings_ch = Channel.of('This', 'is', 'a', 'channel')

// Inspect a channels contents with the operator .view()
strings_ch.view()