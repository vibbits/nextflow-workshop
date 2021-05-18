#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

// Similar to DSL1, the input data is defined in the beginning.
params.reads = "$launchDir/../../data/*{1,2}.fq.gz"
params.outdir = "$launchDir/results"
params.threads = 2
params.slidingwindow = "SLIDINGWINDOW:4:15"
params.avgqual = "AVGQUAL:30"


println """\
      LIST OF PARAMETERS
================================
            GENERAL
Reads            : $params.reads
Output-folder    : $params.outdir/

          TRIMMOMATIC
Threads          : $params.threads
Sliding window   : $params.slidingwindow
Avg quality      : $params.avgqual
"""

// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

include { fastqc as fastqc_raw; fastqc as fastqc_trim } from "${launchDir}/../../modules/fastqc" 
include { trimmomatic } from "${launchDir}/../../modules/trimmomatic"

// Running a workflow with the defined processes here.  
workflow {
	fastqc_raw(read_pairs_ch) 
  trimmomatic(read_pairs_ch)
  fastqc_trim(trimmomatic.out.trim_fq)
}
