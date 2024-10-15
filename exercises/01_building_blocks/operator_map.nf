#!/usr/bin/env nextflow

Channel
    .of( 1, 2, 3, 4, 5 )
    .map { number -> number * number }
    .view()

// older notation you might encounter:
// Channel
//     .of( 1, 2, 3, 4, 5 )
//     .map { it * it }
//     .view()
