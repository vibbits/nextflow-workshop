#!/usr/bin/env nextflow

process python {
    
    script:
    """
    #!/usr/bin/env python3

    firstWord = 'hello'
    secondWord = 'folks'
    print(f'{firstWord} {secondWord}')
    """
}

workflow {
    python()
}

// result visible in work directory in .command.out

