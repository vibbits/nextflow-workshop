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
Output-folder    : ${params.outdir}/

          TRIMMOMATIC
Threads          : ${params.threads}
Sliding window   : ${params.slidingwindow}
Avg quality      : ${params.avgqual}
"""

// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

include { fastqc as fastqc_raw; fastqc as fastqc_trim } from "${projectDir}/../../modules/fastqc" 
include { trimmomatic } from "${projectDir}/../../modules/trimmomatic"

// Running a workflow with the defined processes here.  
workflow {
  read_pairs_ch.view()
  fastqc_raw(read_pairs_ch) 
  trimmomatic(read_pairs_ch)
  fastqc_trim(trimmomatic.out.trim_fq)
}
