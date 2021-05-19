#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process python {
    
    script:
    """
    #!/usr/bin/python3

    firstWord = 'hello'
    secondWord = 'folks'
    print(f'{firstWord} {secondWord}')
    """
}

workflow {
    python()
}

// result visible in work directory in .command.out

