#!/usr/bin/env nextflow

Channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()
