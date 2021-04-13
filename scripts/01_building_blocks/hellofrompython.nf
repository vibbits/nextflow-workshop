#!/usr/bin/env nextflow
 
process python {
    
    """
    #!/usr/bin/python3

    firstWord = 'hello'
    secondWord = 'folks'
    print(f'{firstWord} {secondWord}')
    """
}

// result visible in work directory in .command.out