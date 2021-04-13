
# Building blocks 
In the first chapter we will elaborate on how Nextflow is designed, its advantages and disadvantages, the basic components, etc. 

## Introduction
Writing pipelines to automate processes is not something new. In the `data/` folder we've written a bash script that downloads the data that we will use throughout this tutorial. These bash scripts are probably one of the oldest forms of pipelines where we concatenate processes. Let's have a look at another example:

```
#!/bin/bash

blastp -query sample.fasta -outfmt 6 \
	| head -n 10 \
	| cut -f 2 \
	| blastdbcmd -entry - > sequences.txt
```

Starting with a shebang line, the `blastp` command is piped through multiple times to eventually result in an output file `sequences.txt`. 


```{tab} Question
What is the downside of similar relatively simple pipelines?   
```

```{tab} Solution
There are a couple of suboptimal things happening here:
- Sequential flow of data, there is no parallelization
- Unclear which version of tools are used
- Will it work on my machine? Are the tools installed, is the data there? What about HPC clusters or Cloud environments?
- What if the pipeline fails somewhere in the middle, we need to restart the pipeline from the beginning

```

---

In response to that, workflow managers such as Nextflow were built, aimed to deal with more complex situations. Nextflow is designed around the idea that Linux has many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations. 


By definition, Nextflow is a reactive workflow framework and a programming Domain Specific Language that eases the writing of data-intensive computational pipelines[[1](https://www.nextflow.io/)]. Nextflow scripting is an extension of the Groovy programming language, which in turn is a super-set of the Java programming language. Groovy can be considered as Python for Java in a way that simplifies the writing of code and is more approachable. 

```{image} ../img/nextflow/java-groovy-nextflow.png
:align: center
```

## Why (not)?
Nextflow is not the only player in the field[[2](https://github.com/pditommaso/awesome-pipeline/)], however there are good reasons to opt for it. 

+ Parallelization: processes are automatically scheduled based on available resources 
+ Scalability: simple scaling from local to HPC-cluster usage
+ Portability: run across different platforms
+ Reproducible: native support for containers, conda environments, and interaction with Git. 
+ Re-usability: with the introduction of modules it becomes (theoretically) simple to re-use processes written in other pipelines
+ Community[[3](https://nf-co.re/)]: even though the community is never a reason why to choose for a tool (functionality is more important), it is still very relevant to know that when you are facing problems, there are people out there ready to help you out. 



Some thoughts or disadvantages from my personal point of view. It takes some time to get used to the syntax of the Groovy language. As flexible as it is, as complex it gets. Often it's difficult to trace down the exact problem of a failure of a pipeline script, especially in the beginning. It's probably not the first thing you should be concerned of if you're doing a one-time analysis. 


<!-- Fast prototyping => Custom DSL that enables tasks composition, simplifies most use cases + general purpose programming language for corner cases Easy parallelisation => declarative reactive programming model based on dataflow paradigm, implicit portable parallelism Decouple components => functional approach a task execution is idempotent, ie cannot modify the state of other tasks + isolate dependencies with containers Portable deployments => executor abstraction layer + deployment configuration from implementation logic)
-->

## Main abstractions
Nextflow consists of four main components: channels, operators, processes and workflows. 
- *Channels*: contain the input of the workflows used by the processes. Channels connect processes/operators with each other. 
- *Operators*: transform the content of channels by applying functions or transformations. Usually operators are applied on channels to get the input of a process in the right format.  
- *Processes*: define the piece of script that is actually being run (e.g. an alignment process with STAR). 
- *Workflows*: call the processes as functions with channels as input arguments, only processes defined in the workflow are run. Workflows were introduced in DSL2.


```{image} ../img/nextflow/nextflow-conceptually.png
:align: center
```

The script [`firstscript.nf`](https://github.com/tmuylder/nextflow-workshop/blob/main/scripts/01_building_blocks/firstscript.nf) is using these three components and gives an idea of how Nextflow scripts are being build. 

```bash
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Creating channels
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
```

````{note}
Besides these main building blocks, we also already highlight the existence of the `params` parameters. In the previous code block we explicitly defined some input values in the channels. However, we can define the input values into a parameter instead, that is passed on to the channel. 

```
// create a parameter 'input_read'
params.input_read = '/path/to/read_1.fq'

// use the input_read parameter as an input for the channel
input_read_ch = fromPath(input_read)
```

Here `params.input_read = '/path/to/read_1.fq'` will create a parameter `input_read` and give it the value `'/path/to/read_1.fq'` which is used as an input for the channel. We will later see that these parameters can then be overwritten on runtime. 
````

<!--
(The workflows can be repesented as graphs where the nodes are the processes and the edges are the channels. The processes are block of code that can be executed such as scripts or programs, while the channels are asynchronous queue able to connect processess among them via input / output.)



(Each process is independent from the other and can be run in parallel depending on the availability of processors or if you are in a cluster environment with a scheduler supported by Nextflow. Note also the implicit parallelisation *.fastq in a channel one channel will split it out over multiple processes simultaneously. No need of making a fors–loop.)

(In the previous example the processes A, B and C can be run in parallel and only at their end the process D is triggered.)

-->



### 1. Channels  
The input of the analysis is stored in a channel, these are generally files like sequencing, reference fasta, annotation files, etc. however the input can be of any kind like numbers, strings, lists, etc. To have a complete overview, we refer to the official documentation[[4](https://www.nextflow.io/docs/latest/channel.html#)]. Here are some examples of how a channel is being created:
```
# Channel consisting of strings
strings_ch = Channel.from('This', 'is', 'a', 'channel')

# Channel consisting of a single file
file_ch = Channel.fromPath('data/sequencefile.fastq')

# Channel consisting of multiple files by using a wildcard *
multfiles_ch = Channel.fromPath('data/*.fastq')
```
These channels can then be used by operators or serve as an input for the processes.


````{tab} Exercise 1.1.1
Create a channel consisting of multiple paired-end files by using wildcard * and options {1,2}.
````

````{tab} Solution 1.1.1
```
paired_ch = Channel.fromFilePairs('data/*{1,2}.fastq')

```
````

---




````{tab} Exercise 1.1.2
Which operators do you need to create a channel from a csv-file (`input.csv`)? Write how you would generate the channel. 
````

````{tab} Solution 1.1.2
The file is imported with `.fromPath()`, followed by the `splitCsv()` operator where we set the header to `True`. The last step will output how the channels are constructed. Each row is transformed into a tuple with the first element as a variable `sampleId`, the second as `forward_read` and the third as `reverse_read`.

```
samples_ch = Channel
                .fromPath('input.csv')
                .splitCsv(header:true)
                .view{ row -> tuple(row.sampleId, file(row.forward_read), file(row.reverse_read)) }
```
````

---


### 2. Operators
Operators are necessary to transform the content of channels in a format that is necessary for usage in the processes. There is a plethora of different operators[[5](https://www.nextflow.io/docs/latest/operator.html?highlight=view#)], however only a handful are used extensively. Here are some examples that you might come accross:
- `collect`: e.g. when using a channel consisting of multiple independent files (e.g. fastq-files) and need to be assembled for a next process. 
```
Channel
    .from( 1, 2, 3, 4 )
    .collect()
    .view()

# outputs
[1,2,3,4]
```

- `mix`: e.g. when assembling items from multiple channels into one channel for a next process (e.g. multiqc)

```
c1 = Channel.from( 1,2,3 )
c2 = Channel.from( 'a','b' )
c3 = Channel.from( 'z' )

c1 .mix(c2,c3)

# possible output
a
1
2
b
3
z
```


````{tab} Exercise 1.2.1
Change the last operator in the script from the previous exercise (where we read in the `input.csv` file). Use a transforming operator to prepare the data for further usage. 
````

````{tab} Solution 1.2.1

```
samples_ch = Channel
                .fromPath('input.csv')
                .splitCsv(header:true)
                .map{ row-> tuple(row.sampleId, file(row.forward_read), file(row.reverse_read)) }
```
````

---

### 3. Processes
Processes are the backbone of the pipeline. They represent each individual subpart of the analysis. In the code-snippet below, you can see that it consists of a couple of blocks: directives, input, output, when-clause and the script itself. 

```
process < name > {

   [ directives ]

   input:
    < process inputs >

   output:
    < process outputs >

   when:
    < condition >

   [script|shell|exec]:
   < user script to be executed >
}
```

Each process is executed independently and isolated from any other process. They communicate via asynchronous FIFO queues, i.e. one process will wait for the output of another and then runs reactively when the channel has contents. 



```{image} ../img/nextflow/asynchronous-FIFO.png
:align: center
```

Here is a couple of examples of processes:

````{tab} FastQC
Quality control process with `fastqc`
```
process fastqc {
  input:
  tuple val(sample), path(reads)

  output:
  path("*_fastqc.{zip,html}") 

  script:
  """
  fastqc ${reads}
  """
}
```
````

````{tab} Salmon
Quantifying in mapping-based mode with `salmon`
```
process salmon_quant {
    input:
    path index 
    tuple val(pair_id), path(reads) 

    output:
    path pair_id 

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}
```
````

````{tab} Trimming & quality filtering reads
Trimming adapters & quality filtering with `trimmomatic`
```
process trimmomatic {
    // directives
    publishDir "$params.outdir/trimmed-reads", mode: 'copy', overwrite: true
    label 'low'
    container 'quay.io/biocontainers/trimmomatic:0.35--6'

    input:
    tuple val(sample), path(reads) 

    output:
    tuple val("${sample}"), path("${sample}*_P.fq"), emit: trim_fq
    tuple val("${sample}"), path("${sample}*_U.fq"), emit: untrim_fq
    
    script:
    """
    trimmomatic PE -threads $params.threads ${reads[0]} ${reads[1]} ${sample}1_P.fq ${sample}1_U.fq ${sample}2_P.fq ${sample}2_U.fq $params.slidingwindow $params.avgqual 
    """
}
```
````

---

## Running our first pipeline
If we want to run a Nextflow script in its most basic form, we will use the following command:
```
nextflow run example.nf
```
In our case, we will replace `example.nf` with `02-basic-consepts/firstscript.nf`. First, inspect the script `02-basic-consepts/firstscript.nf` and notice how the channels are being created, passed on to the process' inputs, processed by the script section and then given to the output. 

When we run this script, the result file will not be present in our folder structure. Where will the output of this script be stored?

Nextflow will generate an output that has a standard lay-out:
```
N E X T F L O W  ~  version 20.07.1
Launching `02-basic-concepts/firstscript.nf` [elegant_curie] - revision: 9f886cc00a
executor >  local (2)
executor >  local (2)
[5e/195314] process > valuesToFile (2) [100%] 2 of 2 ✔
results file: /path/to/work/51/7023ee62af2cb4fdd9ef654265506a/result.txt
results file: /path/to/work/5e/195314955591a705e5af3c3ed0bd5a/result.txt
```
The output consists of:
- Version of nextflow 
- Information regarding the script that has ran with an identifier name
- Hash with process ID, progress and caching information
- Optional output printed to the screen as defined in the script (if present)
The results are stored in the results file as described in the two last lines. By default the results of a process are stored in the `work/` directory in subfolders with names defined by the hashes. 

Besides the output, also a bunch of hidden `.command.*` files are present:

- .exitcode, contains 0 if everything is ok, another value if there was a problem.
- .command.log, contains the log of the command execution. Often is identical to .command.out
- .command.out, contains the standard output of the command execution
- .command.err, contains the standard error of the command execution
- .command.begin, contains what has to be executed before .command.sh
- .command.sh, contains the block of code indicated in the process
- .command.run, contains the code made by nextflow for the execution of .command.sh and contains environmental variables, eventual invocations of linux containers etc

Earlier, we described that Nextflow uses an asynchronous FIFO principle. Let's exemplify this by running the script `02-basic-consepts/fifo.nf` and inspect the order that the channels are being processed. 

```
N E X T F L O W  ~  version 20.07.1
Launching `02-basic-concepts/fifo.nf` [nauseous_mahavira] - revision: a71d904cf6
[-        ] process > whosfirst -
This is job number 6
This is job number 3
This is job number 7
This is job number 8
This is job number 5
This is job number 4
This is job number 1
This is job number 2
This is job number 9
executor >  local (10)
[4b/aff57f] process > whosfirst (10) [100%] 10 of 10
```

---

A script, as part of the process, can be written in any language (bash, Python, Perl, Ruby, etc.). This allows to add self-written scripts in the pipeline. The script can be written in the process itself, or can be present as a script in another folder and is run from the process here. 

```
#!/usr/bin/env nextflow
 
process python {
    
    """
    #!/usr/bin/python3

    firstWord = 'hello'
    secondWord = 'folks'
    print(f'{firstWord} {secondWord}')
    """
}
```

```{note}
Directives will be handled further on in the course, conditionals are not considered here.
```



Exercises:

Using the operator view to inspect the result, print out a message in the output of Nextflow running 

```
The results file can be found in: /path/to/...
The results file can be found in: /path/to/...
```

Solution
```
result_ch.view{ "The results file can be found in: $it"  }
```

Exercise:  

You need to execute a task for each record in one or more CSV files.

Strategy:  

Read the CSV file line-by-line using the `splitCsv` operator, then use the `map` operator to return a tuple with the required field for each line. Finally use the resulting channel as input for the process.

Code:   

Given the file `input.csv` with the following content:  

| sampleId | Read 1                        | Read 2                        |
|----------|-------------------------------|-------------------------------|
| 01       | reads/sample_01_read_1.fq.gz  | reads/sample_01_read_2.fq.gz  |
| 02       | reads/sample_02_read_1.fq.gz  | reads/sample_02_read_2.fq.gz  |
| 03       | reads/sample_03_read_1.fq.gz  | reads/sample_03_read_2.fq.gz  |
| 04       | reads/sample_04_read_1.fq.gz  | reads/sample_04_read_2.fq.gz  |


Solution:
```
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_csv = 'exercises/input.csv'

samples_ch = Channel
                .fromPath(params.input_csv)
                .splitCsv(header:true)
                .map{ row-> tuple(row.sampleId, file(row.forward_read), file(row.reverse_read)) }

process split_csv {
    input:
    tuple val(sampleId), file(read1), file(read2)  

    script:
    """
    echo your_command --sample $sampleId --reads $read1 $read2 
    """
}

workflow{
    samples_ch.view()
    split_csv(samples_ch)
}
``` 