#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fw_primer = "GTGCCAGCMGCCGCGGTAA"

process REVERSE_COMP {
    container 'rackspacedot/python37'
    //container 'julia'
    //container 'r-base'
    //container 'debian'

    input:
    val( primer_seq )

    output:
    stdout

    script:
    //PYTHON3
    """
    #!/usr/bin/env python3

    complement_dict = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'U':'A', 'N':'N', 'Y':'R', 'R':'Y', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'V':'B', 'H':'D', 'D':'H'}
    rev_comp = ''.join([complement_dict[x] if x in complement_dict.keys() else 'X' for x in "${primer_seq}"])[::-1]

    print(rev_comp)
    """
}

workflow {
    REVERSE_COMP( params.fw_primer  )
    REVERSE_COMP.out.view()
}