#!/usr/bin/env nextflow

// General parameters
params.datadir = "${launchDir}/data"
params.outdir = "${launchDir}/results"

// Input parameters
params.reads = "${params.datadir}/*{1,2}.fq.gz"
params.genome = "${params.datadir}/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.gtf = "${params.datadir}/ggal_1_48850000_49020000.bed.gff"

// Trimmomatic
params.slidingwindow = "SLIDINGWINDOW:4:15"
params.avgqual = "AVGQUAL:30"

// Star
params.threads = 2
params.genomeSAindexNbases = 10
params.lengthreads = 98

include { fastqc as fastqc_raw; fastqc as fastqc_trim } from "../../../modules/fastqc" //addParams(OUTPUT: fastqcOutputFolder)
include { trimmomatic } from "../../../modules/trimmomatic"
include { star_idx; star_alignment } from "./star"
include { multiqc } from "../../../modules/multiqc"

// Running a workflow with the defined processes here.  
workflow {
    log.info """\
        LIST OF PARAMETERS
    ================================
                GENERAL
    Data-folder      : ${params.datadir}
    Results-folder   : ${params.outdir}
    ================================
        INPUT & REFERENCES 
    Input-files      : ${params.reads}
    Reference genome : ${params.genome}
    GTF-file         : ${params.gtf}
    ================================
            TRIMMOMATIC
    Sliding window   : ${params.slidingwindow}
    Average quality  : ${params.avgqual}
    ================================
                STAR
    Threads          : ${params.threads}
    Length-reads     : ${params.lengthreads}
    SAindexNbases    : ${params.genomeSAindexNbases}
    ================================
    """

    // Channels are being created. 
    def read_pairs_ch = Channel
            .fromFilePairs(params.reads, checkIfExists:true)

    def genome = Channel.fromPath(params.genome)
    def gtf = Channel.fromPath(params.gtf) 

    // QC on raw reads
    fastqc_raw(read_pairs_ch) 
    
    // Trimming & QC
    trimmomatic(read_pairs_ch)
    fastqc_trim(trimmomatic.out.trim_fq)
    
    // Create index
    star_idx(genome, gtf)

    // Combine channels
    def alignment_input = trimmomatic.out.trim_fq
      .combine(star_idx.out.index)
      .combine(gtf)

    alignment_input.view()
    
    // Mapping 
    star_alignment(alignment_input)
    
    // Multi QC on all results
    def multiqc_input = fastqc_raw.out.fastqc_out
      .mix(fastqc_trim.out.fastqc_out)
      .collect()

    multiqc(multiqc_input)
}


// alternative solution using value channels instead:

include { star_alignment as star_alignment_alt_sol } from "../../../modules/star"
workflow alternative_solution {
    // def read_pairs_ch_alt = Channel
    //   .fromFilePairs(params.reads, checkIfExists:true)
    def read_pairs_ch = Channel
            .fromFilePairs(params.reads, checkIfExists:true)

    def genome_alt = Channel.value(file(params.genome))
    def gtf_alt = Channel.value(file(params.gtf))

    // QC on raw reads
    fastqc_raw(read_pairs_ch) 
        
    // Trimming & QC
    trimmomatic(read_pairs_ch)
    fastqc_trim(trimmomatic.out.trim_fq)

    star_idx(genome_alt, gtf_alt)
    star_alignment_alt_sol(trimmomatic.out.trim_fq, star_idx.out.index, gtf_alt)
}
