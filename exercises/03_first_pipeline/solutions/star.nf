process star_idx {
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    path genome
    path gtf
    
    output:
    path "index_dir/", emit: index

    script:
    """
    mkdir index_dir
    
    STAR --runThreadN ${params.threads} \\
      --runMode genomeGenerate \\
      --genomeDir index_dir/ \\
      --genomeFastaFiles ${genome} \\
      --genomeSAindexNbases ${params.genomeSAindexNbases} \\
      --sjdbGTFfile ${gtf}
    """
}

process star_alignment {
    publishDir "${params.outdir}/mapped-reads/", mode: 'copy', overwrite: true  //, pattern: "*.bam"  
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    // (trim_fq, IDX.out, gtf)
    tuple val(sample), path(reads), path(indexDir), path(gtf) 

    output:
    path("*.bam"), emit: align_bam

    script:
    """
    STAR  \\
        --readFilesIn ${reads} \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile ${gtf} \\
        --outFileNamePrefix ${sample}. \\
        --genomeDir ${indexDir}
    """
}


