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
  def x = Channel.value(1)
  def y = Channel.of('a', 'b', 'c')
  foo(x, y)
  foo.out.view()
}