#!/usr/bin/env nextflow

Channel.from( 1,2,3,4, 5 )
    .collect()
    .view()
