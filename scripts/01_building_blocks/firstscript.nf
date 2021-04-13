#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Some comment

/**
 * Some 
 * longer
 * comment
 */

// Creating a channel
numbers_ch = Channel.from(1,2,3)
strings_ch = Channel.from('a','b')

// Defining the process that is executed
process valuesToFile {

    input: 
    val nums
    val strs
    
    output:
    path 'result.txt', emit: result_ch
    
    """
    echo $nums and $strs > result.txt
    """
}

// Running a workflow with the defined processes  
workflow{
    valuesToFile(numbers_ch, strings_ch)
    valuesToFile.out.result_ch.view()
}

