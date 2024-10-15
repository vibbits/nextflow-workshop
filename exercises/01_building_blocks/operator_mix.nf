#!/usr/bin/env nextflow

def c1 = Channel.of( 1,2,3 )
def c2 = Channel.of( 'a','b' )
def c3 = Channel.of( 'z' )

c1.mix(c2,c3)
  .view()