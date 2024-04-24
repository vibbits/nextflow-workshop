#!/usr/bin/env nextflow

process foo {
  input:
  val x
  val y

  output:
  stdout

  script:
  """
  echo $x and $y
  """
}

workflow {
  x = Channel.of(1)
  y = Channel.of('a', 'b', 'c')
  foo(x, y)
  foo.out.view()
}