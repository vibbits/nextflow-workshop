#!/usr/bin/env nextflow

// Some comment

/**
 * Some 
 * longer
 * comment
 */

// Creating a channel
numbers_ch = Channel.of(1,2,3)
strings_ch = Channel.of('a','b')

// Defining the process that is executed
// Defining the process that is executed
process valuesToFile {
    tag  "$nums,$strs" 	

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
workflow {
    valuesToFile(numbers_ch, strings_ch)
    valuesToFile.out.result_ch.view()
}

