#!/usr/bin/env nextflow


// Create the channels
def strings_ch = Channel.of('This', 'is', 'a', 'channel')

// Inspect a channels contents with the operator .view()
strings_ch.view()