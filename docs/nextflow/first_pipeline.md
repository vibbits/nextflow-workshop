


# Creating our first pipeline
In this chapter we will build a basic RNA-seq pipeline consisting of quality controls, trimming of reads and mapping to a reference genome (excl. counting). We will build the pipeline step by step, starting from quality control with FastQC. The figure below was generated with Nextflow and represents the processes that we will build and the overview of the dataflow from the beginning to the end of the workflow. 


```{image} ../img/nextflow/RNAseq.PNG
:align: center
```


## Quality control with `FastQC` 

The following script can be found and run in `exercises/03_first_pipeline/fastqc.nf`.
```bash
#!/usr/bin/env nextflow

params.reads = "${launchDir}/data/*.fq.gz"

/**
 * Quality control fastq
 */

reads_ch = Channel
    .fromPath( params.reads )
    
process fastqc {

    input:
    path read  
    
    script:
    """
    fastqc ${read}
    """
}

workflow {
    fastqc(reads_ch)
}
```

The first line of our script is always a shebang line, declaring the environment where the OS can find the software (i.e. Nextflow). Generally, the input files and parameters of the processes are first assigned into *parameters* which allows flexibility in the pipeline. Input files are then assigned to channels and they serve as input for the process. 

```{note}
- `$launchDir`: The directory from where the script is launched (replaces `$baseDir` in version >20). 
- There is a great flexibility in the Nextflow (Groovy) language: writing of whitespaces, newlines where channels are created,...
```

Let's first run this script with the following command. If you have `htop` installed, keep an eye on the distribution of the workload and notice how Nextflow parallelises the jobs. 
```
nextflow run exercises/03_first_pipeline/fastqc.nf
```

```{note}
The process in `exercises/03_first_pipeline/fastqc.nf` specifies a container, and the `nextflow.config` file in the same folder activates the use of docker. If this directive was not there or docker was not enabled, you would need to make sure that the tool `fastQC` is installed. Conda is already installed and activated, it allows us to easily install `fastqc` with the following command `conda install -c bioconda fastqc`.
```


In the following steps we will add new features to this script:

````{tab} Exercise 2.1
- Overwrite the parameter `reads` on runtime (when running Nextflow on the command-line) so that it only takes `ggal_gut_1.fq.gz` as an input read. 
- Additionally, FastQC generates a html- and zip-file for each read. Where are these output files located? 
````
````{tab} Solution 2.1
```
nextflow run exercises/03_first_pipeline/fastqc.nf --reads data/ggal_gut_1.fq.gz
```
- The output files are stored in the `work/` directory following the generated hashes. The hash at the beginning of each process reveals where you can find the result of each process. 
````

---


````{tab} Exercise 2.2
Change the the script in order to accept & work with paired-end reads. For this we will need to:
- Adapt something in the reads parameter (`params.reads`)
- Change how the channel is generated 
- Change the `input` declaration in the process (from `path` to a `tuple`).  
````
````{tab} Solution 2.2
The solution is given in `exercises/03_first_pipeline/solutions/2.2_fastqc.nf`. Note that if you run this script, only two processes will be launched, one for each paired-end reads dataset.
````

---
````{tab} Exercise 2.3
Run the script with: 
```
nextflow run exercises/03_first_pipeline/fastqc.nf -bg > log
```
What does the `-bg > log` mean? What would the advantage be?
````
````{tab} Solution 2.3
Run in the background and push output of nextflow to the log file. No need of explicitly using nohup, screen or tmux.
````

--- 

````{tab} Exercise 2.4
Check if the files exist ([`checkIfExists`](https://www.nextflow.io/docs/latest/channel.html)) upon creating the channels and invoke an error by running the nextflow script with wrong reads, e.g. 
```
nextflow run exercises/03_first_pipeline/fastqc.nf --reads wrongfilename
````

````{tab} Solution 2.4
The solution is given in `exercises/03_first_pipeline/solutions/2.4_fastqc.nf`
````

--- 

````{tab} Exercise 2.5
Control where and how the output is stored. Have a look at the directive [`publishDir`](https://www.nextflow.io/docs/latest/process.html?highlight=publishdir#publishdir). Nextflow will only store the files that are defined in the `output` declaration block of the process, therefore we now also need to define the output. Put a copy of the output files in a new folder that contains only these results. 


````
````{tab} Solution 2.5
The solution is given in `exercises/03_first_pipeline/solutions/2.5_fastqc.nf` 

- Without any additional arguments, a hyperlink will be created to the files stored in the `work/` directory, with mode set to copy (`mode: 'copy'`) the files will be made available in the defined directory. 
- If the output is to be used by another process, and the files are being moved, they won't be accessible for the next process and hence you're pipeline will fail complaining about files not being present.

```{warning}
Files are copied into the specified directory in an asynchronous manner, thus they may not be immediately available in the published directory at the end of the process execution. For this reason files published by a process must not be accessed by other downstream processes.
```

````
--- 


The final FastQC script, with some additional comments is provided in `exercises/03_first_pipeline/solutions/fastqc_final.nf`. 


## Quality filtering with `trimmomatic` 

Now we will add the next step in our pipeline, which is **trimming and filtering the low quality reads**. For this process, we will use the tool `trimmomatic`. 

The `fastqc.nf` script was extended with the trimmomatic process and is available in `exercises/03_first_pipeline/trimmomatic.nf`. 
- A number of parameters have been added related to the trimmomatic process
- The process `trimmomatic` with its inputs and outputs and the script has been created
- The `workflow` now also contains the process trimmomatic, called as a function

In the `output` declaration block, we are introducing a new option: `emit`. Defining a process output with `emit` allows us to use it as a named channel in the external scope.  

---

At this point we're interested in the result of the `trimmomatic` process. Hence, we want to verify the quality of the reads with another `fastqc` process. Re-run `fastqc` on the filtered read sequences by adding it in the workflow of `trimmomatic.nf`. Use the parameter `-resume` to restart the pipeline from where it stopped the last time. 
  
Hmm, error? `Process fastqc has been already used -- If you need to reuse the same component include it with a different name or include in a different workflow context`. It means that processes can only be used once in a workflow. This means that we need to come up with a smarter solution (see below). 

## Modules
Until now, we have written the processes and the workflow in the same file. However, if we want to be truly modular, we can write a library of modules and import a specific component from that library. A module can contain the definition of a function, process and workflow definitions. 

The figure below gives an overview of how the structure could look like. On the left we have the main Nextflow script (`main.nf`) that defines the parameters, channels and the workflow. It imports the processes from the modules, in this case available in a folder `modules/`. The configuration file `nextflow.config` will be further discussed in the next chapter.  


```{image} ../img/nextflow/overview-folder-structure.png
:align: center
```

A module is generally imported with 
```
include {<process-name>} from './path/to/modules/script.nf'
```
with `<process-name>` the name of the process defined in the `script.nf`. The origin of the module defined by a relative path must start with `./`, alternatively use `projectDir` to use the absolute path. Navigate to the modules folder and find a script called `fastqc.nf`. This script consists of a process and a workflow. This module can be imported into our pipeline script (main workflow) like this:

```
include {fastqc} from './modules/fastqc.nf'
```

This doesn't overcome the problem that we can only use a process once. However, when including a module component it’s possible to specify a name alias. This allows the inclusion and the invocation of the same component multiple times in your script using different names. For example:
``` 
include { fastqc as fastqc_raw; fastqc as fastqc_trim } from "${projectDir}/modules/fastqc"
```
Now we're ready to use a process, defined in a module, multiple times in a workflow. 

Investigate & run the script `exercises/03_first_pipeline/modules.nf` which contains the following code snippet
```
...
include { fastqc as fastqc_raw; fastqc as fastqc_trim } from "${projectDir}/../../modules/fastqc" 
include { trimmomatic } from "${projectDir}/../../modules/trimmomatic"

// Running a workflow with the defined processes here.  
workflow {
  read_pairs_ch.view()
  fastqc_raw(read_pairs_ch) 
  trimmomatic(read_pairs_ch)
  fastqc_trim(trimmomatic.out.trim_fq)
}
```

Similarly as described above, we can extend this pipeline and map our trimmed reads on a reference genome. First, we'll have to create an index for our genome and afterwards we can map our reads onto it. These modules are called from the main script `RNAseq.nf`. 


````{tab} Exercise 2.6
In the folder `modules/` find the script `star.nf` which contains two processes: `star_index` and `star_alignment`. Complete the script `RNAseq.nf` so it includes these processes and hence the pipeline is extended with an indexing and alignment step. The parameters used in the modules are already defined for you. 
````

````{tab} Solution 2.6
Solution in `exercises/03_first_pipeline/solutions/2.6_RNAseq.nf`. The following lines were added. 
```
genome = Channel.fromPath(params.genome)
gtf = Channel.fromPath(params.gtf)

include { star_idx; star_alignment } from "${projectDir}/../../modules/star"

workflow {
  ...
  star_idx(genome, gtf)
  star_alignment(trimmomatic.out.trim_fq, star_idx.out.index, gtf)
}

```` 

---

````{tab} Exercise 2.7
In the folder `modules/` find the script `multiqc.nf`. Import the process in the main script so we can call it as a function in the workflow. This process expects all of the zipped and html files from the fastqc processes (raw & trimmed) as one input. Thus it is necessary to use the operators `.mix()` and `.collect()` on the outputs of `fastqc_raw` and `fastqc_trim` to generate one channel with all the files.  
````

````{tab} Solution 2.7
Solution in `exercises/03_first_pipeline/solutions/2.7_RNAseq.nf`. The following lines were added. 
```
include { multiqc } from "${projectDir}/../../modules/multiqc" 

workflow {
  ...
  multiqc_input = fastqc_raw.out.fastqc_out
    .mix(fastqc_trim.out.fastqc_out)
    .collect()

  multiqc(multiqc_input)
} 
```
```` 

---

You might have noticed that the star_alignment process was only executed once in exercise 2.6 and 2.7, while we expect the process to be executed twice (we have 2 samples). This is due to the way we have defined the input for the star_alignment process.

```bash
process star_alignment {
    publishDir "${params.outdir}/mapped-reads/", mode: 'copy', overwrite: true  //, pattern: "*.bam"  
    label 'high'
    container "quay.io/biocontainers/star:2.6.1d--0"

    input:
    tuple val(sample), path(reads) 
    path indexDir
    path gtf

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
```
As you can see, we have defined 3 separate input channels for our process.


When two or more channels are declared as process inputs, the process waits until there is a complete input configuration, i.e. until it receives a value from each input channel. When this condition is satisfied, the process consumes a value from each channel and launches a new task, repeating this logic until one or more channels are empty.
More information can be found in the [documentation](https://www.nextflow.io/docs/latest/process.html#multiple-input-channels)

Because we have more than 1 sample in the first input channel, but only 1 entry for both the second (indexDir) and third (gtf) channel, the process will only be executed once.

---

````{tab} Exercise 2.8
Find a way to restructure the input channel for the `star_alignment` process so it will correctly be exectuted for each sample instead of just once. 

- Use channel operators to combine the multiple input channels 
- Don't forget to change the input declaration in the process as well

````

````{tab} Solution 2.8
Solution in `exercises/03_first_pipeline/solutions/2.8_RNAseq.nf`. The following lines were added. 

```
workflow {
  ...
  // Combine channels
  alignment_input = trimmomatic.out.trim_fq
    .combine(star_idx.out.index)
    .combine(gtf)

  alignment_input.view()

  // Mapping 
  star_alignment(alignment_input)
}
```
The following adjustments were made to the input declaration block of the `star.nf` module.
```
process star_alignment {
    ...
    input:
    // (trim_fq, IDX.out, gtf)
    tuple val(sample), path(reads), path(indexDir), path(gtf) 

    ...
}


```
```` 


---

This pipeline is still subject to optimizations which will be further elaborated in the next chapter. 


## Subworkflows

The workflow keyword allows the definition of **sub-workflow** components that enclose the invocation of one or more processes and operators. Here we have created a sub-workflow for a hypothetical `hisat` aligner. 

```
workflow hisat {
  hisat_index(arg1)
  hisat_alignment(arg1, arg2)
}
```

These sub-workflows allow us to use this workflow from within another workflow. The workflow that does not cary any name is considered to be the main workflow and will be executed implicitly. This is thus the entry point of the pipeline, however alternatively we can overwrite it by using the `-entry` parameter. The following code snippet defines two sub-workflows and one main workflow. If we would only be interested in the star alignment workflow, then we would use `nextflow run pipeline.nf -entry star`. 

```
workflow star {
  take:
  arg1
  arg2
  arg3

  main:
  star_index(arg1, arg2)
  star_alignment(arg1, arg2, arg3)
}

workflow hisat2 {
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

```{note}
The `take:` declaration block defines the input channels of the sub-workflow, `main:` is the declaration block that contains the processes (functions) and is required in order to separate the inputs from the workflow body. These options are useful when the pipeline is growing with multiple entry-levels to keep a tidy overview. 
```




## Extra exercises

````{tab} Extra exercise 1
Extend the workflow pipeline with a final note printed on completion of the workflow. Read more about workflow introspection [here](https://www.nextflow.io/docs/latest/metadata.html). 
````
````{tab} Solution 1
The solution is given in `exercises/03_first_pipeline/solutions/ex.1_RNAseq.nf`
````

---

````{tab} Extra exercise 2
Adapt the `exercises/03_first_pipeline/solutions/ex.1_RNAseq.nf` script so it uses Salmon as an aligner and quantifier. In our temporary solution the alignment with Star has been replaced with Salmon, it would be better to create a subworkflow so you can choose upon `-entry` to work with Star or Salmon.  
````
````{tab} Solution 2
The solution is given in `exercises/03_first_pipeline/solutions/ex.2_RNAseq.nf`. 
````



---

````{tab} Extra exercise 3
Write a Nextflow script for a tool that you use in your research. Use the same approach with parameters, channels, process in a module, and a workflow. 
````
````{tab} Solution 3
If you are stuck, don't hesitate to ask for help! 
````
