#!/usr/bin/env nextflow

Channel
    .of( 1, 2, 3, 4, 5 )
    .map { it * it }
    .subscribe onNext: { println it }, onComplete: { println 'Done' }

// using .view() to vieuw the channel:
// Channel
//     .of( 1, 2, 3, 4, 5 )
//     .map { it * it }
//     .view()


// same result with different notation:
// Channel
//     .of( 1, 2, 3, 4, 5 )
//     .map { it -> it * it }
//     .view()

// same result with different notation:
// Channel
//     .of( 1, 2, 3, 4, 5 )
//     .map { bar -> bar * bar }
//     .view()