#!/usr/bin/env nextflow

Channel
    .of( 1, 2, 3, 4, 5 )
    .map { it * it }
    .subscribe onNext: { println it }, onComplete: { println 'Done' }

// Channel
//     .of( 1, 2, 3, 4, 5 )
//     .map { it * it }
//     .view()