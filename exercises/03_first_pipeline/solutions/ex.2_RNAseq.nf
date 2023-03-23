#!/usr/bin/env nextflow



// General parameters
params.datadir = "${launchDir}/data"
params.outdir = "${launchDir}/results"

// Input parameters
params.reads = "${params.datadir}/*{1,2}.fq.gz"
params.transcriptome = "${params.datadir}/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"

// Trimmomatic
params.slidingwindow = "SLIDINGWINDOW:4:15"
params.avgqual = "AVGQUAL:30"

// Salmon
params.threads = 2


log.info """\
      LIST OF PARAMETERS
================================
            GENERAL
Data-folder      : ${params.datadir}
Results-folder   : ${params.outdir}
================================
      INPUT & REFERENCES 
Input-files      : ${params.reads}
Reference genome : ${params.transcriptome}
================================
          TRIMMOMATIC
Sliding window   : ${params.slidingwindow}
Average quality  : ${params.avgqual}
================================
             SALMON
Threads          : ${params.threads}
================================
"""


// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

transcriptome = Channel.fromPath(params.transcriptome)

include { fastqc as fastqc_raw; fastqc as fastqc_trim } from "${projectDir}/../../../modules/fastqc" //addParams(OUTPUT: fastqcOutputFolder)
include { trimmomatic } from "${projectDir}/../../../modules/trimmomatic"
include { salmon_idx ; salmon_quant } from "${projectDir}/../../../modules/salmon"
include { multiqc } from "${projectDir}/../../../modules/multiqc" 

// Running a workflow with the defined processes here.  
workflow {
  // QC on raw reads
  fastqc_raw(read_pairs_ch) 
	
  // Trimming & QC
  trimmomatic(read_pairs_ch)
  fastqc_trim(trimmomatic.out.trim_fq)
	
  // Mapping salmon
  salmon_idx(transcriptome)
  salmon_quant(salmon_idx.out, read_pairs_ch)
  
  // Multi QC on all results
  multiqc_input = fastqc_raw.out.fastqc_out
    .mix(fastqc_trim.out.fastqc_out)
    .collect()

  multiqc(multiqc_input)
}

workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}