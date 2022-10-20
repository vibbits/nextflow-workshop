#!/usr/bin/env nextflow


// Similar to DSL1, the input data is defined in the beginning.
params.reads = "${launchDir}/data/*{1,2}.fq.gz"
params.outdir = "${launchDir}/results"
params.threads = 2
params.slidingwindow = "SLIDINGWINDOW:4:15"
params.avgqual = "AVGQUAL:30"


log.info """\
      LIST OF PARAMETERS
================================
            GENERAL
Reads            : ${params.reads}
Output-folder    : ${params.outdir}

          TRIMMOMATIC
Threads          : ${params.threads}
Sliding window   : ${params.slidingwindow}
Avg quality      : ${params.avgqual}
"""

// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

// Definition of a process, notice the absence of the 'from channel'.
// A process being defined, does not mean it's invoked (see workflow)
process fastqc {
  publishDir "${params.outdir}/quality-control-${sample}/", mode: 'copy', overwrite: true
  container 'quay.io/biocontainers/fastqc:0.11.9--0'
  
  input:
  tuple val(sample), path(reads)

  script:
  """
  fastqc ${reads}
  """
}

// Process trimmomatic
process trimmomatic {
  publishDir "${params.outdir}/trimmed-reads-${sample}", mode: 'copy'
  container 'quay.io/biocontainers/trimmomatic:0.35--6'

  // Same input as fastqc on raw reads, comes from the same channel. 
  input:
  tuple val(sample), path(reads) 

  output:
  tuple val("${sample}"), path("${sample}*_P.fq"), emit: paired_fq
  tuple val("${sample}"), path("${sample}*_U.fq"), emit: unpaired_fq

  script:
  """
  trimmomatic PE -threads ${params.threads} ${reads[0]} ${reads[1]} ${sample}1_P.fq ${sample}1_U.fq ${sample}2_P.fq ${sample}2_U.fq ${params.slidingwindow} ${params.avgqual} 
  """
}

// Running a workflow with the defined processes here.  
workflow {
	fastqc(read_pairs_ch) 
  trimmomatic(read_pairs_ch)
  fastqc(trimmomatic.out.paired_fq) // This will raise an error. Do you remember why?
}
