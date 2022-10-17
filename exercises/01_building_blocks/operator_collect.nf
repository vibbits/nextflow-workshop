#!/usr/bin/env nextflow

Channel
    .from( 1, 2, 3, 4 )
    .collect()
    .view()
