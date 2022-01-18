#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// set default input parameters (these can be altered by calling their flag on the command line, e.g., nextflow run main.nf --reads 'data2/*_R{1,2}.fastq')
params.reads = 'data1/*_R{1,2}.fastq'
params.fw_primer = "GTGCCAGCMGCCGCGGTAA" 
params.rv_primer = "GGACTACHVHHHTWTCTAAT"


// Set a header made using https://patorjk.com/software/taag (but be sure to escape characters such as dollar signs and backslashes, e.g., '$'=> '\$' and '\' =>'\\')
log.info """
    ==============================================================================================

                                            \$\$\\                     \$\$\\ \$\$\\                     
                                            \\__|                    \$\$ |\\__|                    
    \$\$\$\$\$\$\\\$\$\$\$\\  \$\$\\   \$\$\\        \$\$\$\$\$\$\\  \$\$\\  \$\$\$\$\$\$\\   \$\$\$\$\$\$\\  \$\$ |\$\$\\ \$\$\$\$\$\$\$\\   \$\$\$\$\$\$\\  
    \$\$  _\$\$  _\$\$\\ \$\$ |  \$\$ |      \$\$  __\$\$\\ \$\$ |\$\$  __\$\$\\ \$\$  __\$\$\\ \$\$ |\$\$ |\$\$  __\$\$\\ \$\$  __\$\$\\ 
    \$\$ / \$\$ / \$\$ |\$\$ |  \$\$ |      \$\$ /  \$\$ |\$\$ |\$\$ /  \$\$ |\$\$\$\$\$\$\$\$ |\$\$ |\$\$ |\$\$ |  \$\$ |\$\$\$\$\$\$\$\$ |
    \$\$ | \$\$ | \$\$ |\$\$ |  \$\$ |      \$\$ |  \$\$ |\$\$ |\$\$ |  \$\$ |\$\$   ____|\$\$ |\$\$ |\$\$ |  \$\$ |\$\$   ____|
    \$\$ | \$\$ | \$\$ |\\\$\$\$\$\$\$\$ |      \$\$\$\$\$\$\$  |\$\$ |\$\$\$\$\$\$\$  |\\\$\$\$\$\$\$\$\\ \$\$ |\$\$ |\$\$ |  \$\$ |\\\$\$\$\$\$\$\$\\ 
    \\__| \\__| \\__| \\____\$\$ |      \$\$  ____/ \\__|\$\$  ____/  \\_______|\\__|\\__|\\__|  \\__| \\_______|
                  \$\$\\   \$\$ |      \$\$ |          \$\$ |                                            
                  \\\$\$\$\$\$\$  |      \$\$ |          \$\$ |                                            
                   \\______/       \\__|          \\__|                                                  
    
    ==============================================================================================

    INPUT PARAMETERS:
        - reads : ${params.reads}
        - forward primer sequence : ${params.fw_primer}
        - reverse primer sequence : ${params.rv_primer}

    ==============================================================================================
    """.stripIndent()



// set input data
Channel.fromFilePairs( params.reads , checkIfExists:true)
    .set{ pe_reads_ch }


// use groovy logic to calculate the reverse complements of the primer sequences
def complements = [ 
            A:'T', T:'A', U:'A', G:'C',
            C:'G', Y:'R', R:'Y', S:'S', 
            W:'W', K:'M', M:'K', B:'V',
            D:'H', H:'D', V:'B', N:'N'
            ]
def fw_rc_primer = params.fw_primer.collect{ base -> return complements[base] ?: 'X' }.reverse().join()
def rv_rc_primer = params.rv_primer.collect{ base -> return complements[base] ?: 'X' }.reverse().join()


// include processes and subworkflows to make them available for use in this script 
include { check_QC as check_QC_raw; check_QC as check_QC_trimmed } from "${launchDir}/modules/QC" 
include { CUTADAPT } from "${launchDir}/modules/trimming"
include { DADA2 } from "${launchDir}/modules/reads2counts"


//set the path to the script to run in the DADA2 process (you can also make a folder 'bin' and put this script in there so it will automatically be added to nextflow's path)
params.script1 = "${launchDir}/reads2counts.r"


workflow {
    check_QC_raw( "raw", pe_reads_ch ) //pass the 'state' and the raw reads to the QC subworkflow

    CUTADAPT( pe_reads_ch, params.fw_primer, params.rv_primer, fw_rc_primer, rv_rc_primer) //pass the raw reads and the primer sequences to the cutadapt process

    check_QC_trimmed( "trimmed", CUTADAPT.out) //pass the 'state' and the trimmed reads to the QC subworkflow
    
    DADA2( CUTADAPT.out.collect{ x -> x[1] }, params.script1) //pass the paths to the reads and the R script to the DADA2 process
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Time to complete workflow execution: $workflow.duration"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: $workflow.errorMessage"
}