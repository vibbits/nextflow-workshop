#!/usr/bin/env nextflow

// set default input parameters (these can be altered by calling their flag on the command line, e.g., nextflow run main.nf --reads 'data2/*_R{1,2}.fastq')
params.reads = 'data1/*_R{1,2}.fastq'
params.outdir = "${launchDir}/output"
params.fw_primer = "GTGCCAGCMGCCGCGGTAA"
params.fw_primer_rev_comp = "TTACCGCGGCKGCTGGCAC"
params.rv_primer = "GGACTACHVHHHTWTCTAAT"
params.rv_primer_rev_comp = "ATTAGAWADDDBDGTAGTCC"

//set the path to the script to run in the DADA2 process (you can also make a folder 'bin' and put this script in there so it will automatically be added to nextflow's path)
params.script1 = "${projectDir}/reads2counts.r"

// include processes and subworkflows to make them available for use in this script 
include { check_QC as check_QC_raw; check_QC as check_QC_trimmed } from "./modules/QC" 
include { CUTADAPT } from "./modules/trimming"
include { DADA2 } from "./modules/reads2counts"


workflow {
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
        - output directory : ${params.outdir}
        - forward primer sequence : ${params.fw_primer}
        - forward primer reverse complement sequence : ${params.fw_primer_rev_comp}
        - reverse primer sequence : ${params.rv_primer}
        - reverse primer reverse complement sequence : ${params.rv_primer_rev_comp}

    ==============================================================================================
    """.stripIndent()

    // set input data
    def pe_reads_ch = Channel
        .fromFilePairs(params.reads , checkIfExists:true)

    //pass the 'step' and the raw reads to the QC subworkflow
    check_QC_raw("raw", pe_reads_ch)
    
    // the "raw" notation creates a value channel. This is equivalent to the following lines
    // step1 = channel.value("raw")
    // check_QC_raw(step1, pe_reads_ch)

    //pass the raw reads and the primer sequences to the cutadapt process
    CUTADAPT(pe_reads_ch)

    //pass the 'step' and the trimmed reads to the QC subworkflow
    check_QC_trimmed("trimmed", CUTADAPT.out)
    
    //pass the paths to the reads to the DADA2 process
    def dada2_input = CUTADAPT.out
        .map{_sample, reads -> reads}
        .collect()
    

    // you could also add the closure to the collect operator to do this in one step
    // dada2_input = CUTADAPT.out
    //     .collect{x -> x[1]}

    DADA2(dada2_input)

}

workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
