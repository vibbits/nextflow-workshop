


# Creating our first pipeline
In this chapter we will build a basic RNA-seq pipeline consisting of quality controls, trimming of reads and mapping to a reference genome (excl. counting). We will build the pipeline step by step, starting from quality control with FastQC. We will start with DSL1 for this process. After exploring Nextflow's flexibility on the quality control process, we will switch to the newer DSL2 language and extend our pipeline script with the other processes. 


```{image} ../img/nextflow/RNAseq.PNG
:align: center
```

## Quality control with `FastQC` (DSL1)

The following script can be found and run in `03-first-pipeline/fastqc_1.nf`:
```bash
#!/usr/bin/env nextflow

params.reads = "$launchDir/data/*.fq.gz"

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromPath( params.reads )
    .view()
    
process fastqc_raw_reads {

    input:
    file read from reads_ch 
    
    script:
    """
    fastqc ${read}
    """
}
```

The first line of our script is always a shebang line, declaring the environment where the OS can find the software (i.e. Nextflow). Generally, the input files and parameters of the processes are first assigned into *parameters* which allows flexibility in the pipeline. Input files are then assigned to channels and they serve as input for the process. 

```{note}
- `$launchDir`: The directory where the main script is located (replaces `$baseDir` in version >20). 
- There is a great flexibility in the Nextflow (Groovy) language: writing of whitespaces, newlines where channels are created, assigning channel values to a variable or using `.set{}`, etc.) 
```

Let's first run this script with the following command. If you have `htop` installed, keep an eye on the distribution of the workload and notice how Nextflow implicitly parallelises the jobs. 
```
nextflow run 03-first-pipeline/fastqc_1.nf
```

In the following steps we will add new features to this script:
1. Run script with the following line: `nextflow run 03-first-pipeline/fastqc_1.nf -bg > log`. What does the `-bg > log` mean? What would the advantage be? Additionally, FastQC generates an html- and zip-file for each read. Where are the output files located?    

2. Adapt file for handling read pairs. Which parameter would you use on runtime to overwrite the inputfiles? 

3. Print parameters using `println` & check if the files exist when creating the channels. Additionally, invoke an error by running the nextflow script with wrong reads, e.g. `nextflow run 03-first-pipeline/fastqc_3.nf --reads wrongfilename`.
```{hint}
[`checkIfExists`](https://www.nextflow.io/docs/latest/channel.html?highlight=fromfilepairs)
``` 

4. Create a directory where the files can be stored. 

```{hint}
Hint: [`publishDir`](https://www.nextflow.io/docs/latest/process.html?highlight=publishdir#publishdir). 
```

```{note}
  - Without any additional arguments, a hyperlink will be created to the files stored in the `work/` directory, with mode set to copy (`mode: 'copy'`) the files will be made available in the defined directory. 
  - Can you guess what might happen if we set the mode to move? (`mode: 'move'`)  
  - Can you guess what might happen if we set `pattern` to `'*.zip'`?
```

<!--
(ANSWER 1.1: run in the background and push output of nextflow to the log file. No need of explicitly using nohup, screen or tmux.)

(ANSWER 1.2: the hash at the beginning of each process reveals where you can find the result of each process.)  
-->
<!--
(ANSWER 2: Result: [`fastqc_2.nf`]
nextflow run --reads data/WT_lib1_R1.fq.gz 03-first-pipeline/fastqc_1.nf .)  
--> 

<!--
(ANSWER 3: Result: [`fastqc_3.nf`])  
--> 


<!--
(ANSWER 4: Result: [`fastqc_4.nf`]. If the output is to be used by another process, and the files are being moved, they won't be accessible for the next process and hence you're pipeline will fail complaining about files not being present.) 

(ANSWER: Only zip files will be created to the output. Patterns will specify what will be outputted. ie zip and html. So we can control which files are being outputted in the publishdir.) 
--> 

**Warning**: Files are copied into the specified directory in an asynchronous manner, thus they may not be immediately available in the published directory at the end of the process execution. For this reason files published by a process must not be accessed by other downstream processes.

## Moving towards DSL2
Nextflow recently went through a big make-over. The premise of the next version, using DSL2, is to make the pipelines more modular and simplify the writing of complex data analysis pipelines. 

Here is a list of the major changes: 
- Following the shebang line, the nf-script wil start with the following line: `nextflow.enable.dsl=2` (not to be mistaken with *preview*).
- When using DSL1 each channel could only be consumed once, this is ommited in DSL2. Once created, a channel can be consumed indefinitely. 
- A process on the other hand can still only be used once in DSL2 
- A new term is introduced: `workflow`. In the workflow, the processes are called as functions with input arguments being the channels. 
- Regarding the processes, the new DSL separates the definition of a process from its invocation. This means that in DSL1 the process was defined and also run when the script was invoked, however in DSL2, the definition of a process does not necessarily mean that it will be run. 
- Moreover, within processes there are no more references to channels (i.e. `from` and `into`). The channels are passed as inputs to the processes which are defined and invoked in the `workflow`. 


## Quality control with `FastQC` (DSL2)

Let's have a look (`nextflow run 03-first-pipeline/dsl2-fastqc.nf`):

```
#!/usr/bin/env nextflow

// This is needed for activating the new DLS2
nextflow.enable.dsl=2

// Similar to DSL1, the input data is defined in the beginning.
params.reads = "$launchDir/data/*{1,2}.fq.gz"
params.outdir = "$launchDir/results"

println """\
      LIST OF PARAMETERS
================================
Reads            : $params.reads
Output-folder    : $params.outdir/
"""

// Also channels are being created. 
read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists:true)

// Definition of a process, notice the absence of the 'from channel'.
// A process being defined, does not mean it's invoked (see workflow)
process fastqc {
  publishDir "$params.outdir/quality-control-$sample/", mode: 'copy', overwrite: true
    
  input:
  tuple val(sample), file(reads)

  script:
  """
  mkdir -p $params.outdir/quality-control-$sample
  fastqc --outdir $params.outdir/quality-control-$sample ${reads}
  """
}

// Running a workflow with the defined processes here.  
workflow {
	read_pairs_ch.view()
	fastqc(read_pairs_ch) 
}
```

## Trimming with `trimmomatic` (DSL2)

Now we will add the next step in our pipeline, which is **trimming and filtering the low quality reads**. For this process, we will use the tool `trimmomatic`. 

The solution is available in `03-first-pipeline/dsl2-trimming.nf`. Here we're introducing a new option: `emit`. Defining a process output with `emit` allows us to use it as a channel in the external scope.  

At this point we're interested in the result of the `trimmomatic` process. Hence, we want to verify the quality of the reads with another `fastqc` process. Re-run `fastqc` on the filtered read sequences by adding it in the workflow of `03-first-pipeline/dsl2-trimming.nf`. Use the parameter `-resume` to restart the pipeline from where it stopped the last time. 
  
Hmm, error? `Process fastqc has been already used -- If you need to reuse the same component include it with a different name or include in a different workflow context`. It means that processes can only be used once in a workflow. This means that we need to come up with a smarter solution (see below). 

## Subworkflows and modules
The workflow keyword allows the definition of **sub-workflow** components that enclose the invocation of one or more processes and operators. It also allows you to use this workflow from within another workflow. The workflow that does not cary any name is considered to be the main workflow and will be executed implicitly. 

```
workflow star{
  take:
  arg1
  arg2
  arg3

  main:
  star_index(arg1, arg2)
  star_alignment(arg1, arg2, arg3)
}

workflow hisat2{
  take:
  arg1
  arg2

  main:
  hisat_index(arg1)
  hisat_alignment(arg1, arg2)
}

workflow {
  star(arg1, arg2, arg3)
  hisat2(arg1, arg2)
}

```


However, if we want to be truly modular, we can write a library of modules and import a specific component from that library. A module can contain the definition of a function, process and workflow definitions. Navigate to the modules folder and find a script called `fastqc.nf`. This script consists of a process and a workflow. This module can be imported into our pipeline script (main workflow) like this:
```
include {QC} from './modules/fastqc.nf'
```
This line is quite specific. The workflow is defined within the curly brackets, the origin of the module defined by a relative path must start with `./`. 

When including a module component it’s possible to specify a name alias. This allows the inclusion and the invocation of the same component multiple times in your script using different names. For example:
``` 
include { QC as fastqc_raw; QC as fastqc_trim } from "${launchDir}/modules/fastqc"
```
Now we're ready to use a process, defined in a module, multiple times in a workflow. 

Run `03-first-pipeline/dsl2-subworkflow.nf` which contains the `trimmomatic` process internally and imports the `fastqc` process from the modules library `module/fastqc.nf`. 

## RNAseq pipeline
Similarly as described above, we can extend this pipeline and map our trimmed reads on a reference genome. First, we'll have to index our reads and afterwards we can map our reads. In the folder `modules/` find the script `star.nf` which contains two processes: `star_index` and `star_alignment`. 
These modules are called from the main script `03-first-pipeline/dsl2-RNAseq.nf`. 

This pipeline is still subject to optimizations which will be elaborated in the next section. 

