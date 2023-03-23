# Executing pipelines


## Executing our first pipeline
If we want to run a Nextflow script in its most basic form, we will use the following command:
```
nextflow run <pipeline-name.nf>
```

with `<pipeline-name.nf>` the name of our pipeline, e.g. `exercises/02_run_first_script/firstscript.nf`. Inspect the script `firstscript.nf` again and notice how the channels and process are being created, how the workflow calls the process as a function with the channels as input arguments, how they are passed on as the processes' inputs, to the script section and then given to the output. 

```
#!/usr/bin/env nextflow

// Creating a channel
numbers_ch = Channel.of(1,2,3)
strings_ch = Channel.of('a','b')

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
workflow {
    valuesToFile(numbers_ch, strings_ch)
    valuesToFile.out.result_ch.view()
}
```


Nextflow will generate an output that has a standard lay-out:
```
N E X T F L O W  ~  version 22.04.5
Launching `exercises/02_run_first_script/firstscript.nf` [distracted_almeida] DSL2 - revision: 1a87b5fe26
executor >  local (2)
[eb/9af3b0] process > valuesToFile (2) [100%] 2 of 2 ✔
/home/training/git/nextflow-workshop/work/c8/b5f6c2d2a5932f77d5bc53320b8a5d/result.txt
/home/training/git/nextflow-workshop/work/eb/9af3b0384ef96c011b4da69e86fca7/result.txt
```

The output consists of:
- Version of nextflow 
- Information regarding the script that has ran with an identifier name
- Hash with process ID, progress and caching information
- Optional output printed to the screen as defined in the script (if present)

```{admonition} Question
:class: hint
When we run this script, the result file will not be present in our folder structure. Where will the output of this script be stored?
```

The results are stored in the results file as described in the two last lines. By default the results of a process are stored in the `work/` directory in subfolders with names defined by the hashes. Besides the output that we generated, also a bunch of hidden `.command.*` files are present in the hashed `work` folders:

```
|- work/
|   |
|   |- c8
|   |   |
|   |   |- b5f6c2d2a5932f77d5bc53320b8a5d
|   |   |   |
|   |   |   |- .command.begin
|   |   |   |- .command.err
|   |   |   |- .command.log
|   |   |   |- ...
|   |   
|   |- eb
|   |   |
|   |   |- 9af3b0384ef96c011b4da69e86fca7
|   |   |   |
|   |   |   | -...
... 
```

```{tab} .command.log
`.command.log`, contains the log of the command execution. Often is identical to .command.out
```
```{tab} .command.out
`.command.out`, contains the standard output of the command execution
```
```{tab} .command.err
`.command.err`, contains the standard error of the command execution
```
```{tab} .command.begin
`.command.begin`, contains what has to be executed before .command.sh
```
```{tab} .command.sh
`.command.sh`, contains the block of code indicated in the process script block
```
```{tab} .command.run
`.command.run`, contains the code made by nextflow for the execution of .command.sh and contains environmental variables, eventual invocations of linux containers etc
```
```{tab} .exitcode
`.exitcode`, contains the exitcode of the proccess, this is typically 0 if everything is ok, another value if there was a problem.
```

---




## Pipeline parameters vs Nextflow options

There are two types of parameters! 

Pipeline parameters are the parameters used in the pipeline script (e.g. `params.reads`). They are related to the pipeline and can be modified/overwritten on the command-line with a **double dash**: e.g parameter `params.reads` in the `fastqc.nf` script can be set as `--reads` in the command-line. 

There are more ways to set your pipeline parameters, for example in a `nextflow.config` file. This can be useful when there are many parameters to a pipline, or if you want to save the parameters for reuse later. More information about this can be found [here](https://www.nextflow.io/docs/latest/config.html).

Nextflow  options are set in the command-line with a **single dash** and are predefined in Nextflow's language. Here are some examples:
- `-bg` runs the workflow in the background.
- `-resume` resumes the parameter from where it failed last time and uses cached information from the `work/` directory.
- `-with-report` creates a report of how the pipeline ran (performance, memory usages etc.).
- `-work-dir` overwrite the name of the directory where intermediate result files are written.
- ...   

We will discover these options while going through the course materials. 


## Knowing where to find a pipeline and which one to use.
Before thinking of writing our own (plausibly) complex pipeline, we can also think about importing one. Several repositories exist that store Nextflow pipelines (non-exhaustive list):  
    - Some curated nextflow pipelines are available on [awesome-nextflow](https://github.com/nextflow-io/awesome-nextflow).  
    - Pipelines from the [nf-core community](https://nf-co.re/pipelines).  
    - Pipelines from [WorkflowHub](https://workflowhub.eu/) (this is a currently ongoing effort).  
    - VSN-Pipelines for single cell analysis [VSN-Pipelines](https://github.com/vib-singlecell-nf/vsn-pipelines) (Currently not updated)
       

## Import a pipeline 

Imagine that we set our eyes on the [`nextflow-io/rnaseq-nf`](https://github.com/nextflow-io/rnaseq-nf) pipeline. A toy workflow for the analysis of (once again) RNAseq data.

There are different possibilities to pull a publicly available pipeline at a git-based hosting code system (GitHub, GitLab or BitBucket). 
One of them is to pull the pipeline using `nextflow pull`, like so:
```
nextflow pull nextflow-io/rnaseq-nf
```
The latest version of the pipeline is written in DSL2. Imagine that you would like to run the last DSL1 version of the pipeline (v1.2), we can pull this specific version using:
```
nextflow pull nextflow-io/rnaseq-nf -r v1.2
```
Nextflow enables to pull any specific tag, release or commit. To pull the pipeline from (1) a given branch and at a (2) specific git commit, we use the following:
```
nextflow pull nextflow-io/rnaseq-nf -r master
nextflow pull nextflow-io/rnaseq-nf -r 98ffd10a76
```

The workflows will not be cloned in the folder from where we launched these commands. Instead, it is available in the folder `~/.nextflow/assets/`, e.g. for the nextflow-io/rnaseq-nf pipeline in `~/.nextflow/assets/nextflow-io/rnaseq-nf/`. If we would want to have the workflows available (for further editing), we can use `nextflow clone`, similar to how `git` works. 

The `-r` option can also be used directly with `nextflow run` rather than running `nextflow pull` first.


--- 

After importing our pipeline of interest, we can run it on the command-line using the nextflow run `<pipeline-name>` command, with `<pipeline-name>` being the name of the pipeline we just imported. 

```{note}  
When you use `nextflow run` without pulling the pipeline first (`nextflow pull`), Nextflow will check GitHub for a corresponding repository, if one exists it will pull it and run it locally. 

`nextflow run nextflow-io/rnaseq-nf` will result in an error due to uninstalled tools on our VMs. To fix this, simply add the parameter `-with-docker`. We will discover what is happening when we enable this setting later`. 
``` 

---


## Extra exercises

```{warning}
The `nextflow-io/rnaseq-nf` pipeline was recently upgraded to require Nextflow-22.08.2-edge which changes how the `conda` directive works. As we are on an earlier version of Nextflow, you will need to use the version in this git commit `1ca363c8`.
```

````{tab} Extra exercise 1
Run the publicly available pipeline `nextflow-io/rnaseq-nf`. Try to modify the name of the folder where results are stored by using a different parameter on the command-line.
````

````{tab} Solution 1
The directory with the final results: 
```bash
nextflow run nextflow-io/rnaseq-nf --outdir 'myAwesomeResults' -r 1ca363c8 -with-docker
```
or, the directory with temporary files (used for caching): 
```bash
nextflow run nextflow-io/rnaseq-nf -w 'myAwesomeResults' -r 1ca363c8 -with-docker
```
````

--- 

````{tab} Extra exercise 2
Which pipeline parameters are defined, can you modify these in the rnaseq-nf pipeline?
````

````{tab} Solution 2
The `reads`, `transcriptome`, `outdir` and `multiqc` parameters.
````

--- 

````{tab} Extra exercises 3
3.1 How many pipelines are currently available in [nf-core](https://nf-co.re/)? How many are under development, released, and archived?

3.2 Find the pipeline for performing ATAC-seq data analysis in [nf-core](https://nf-co.re/). 
- What is the current/latest version of the pipeline? 
- How many versions are available to download? 
- How many and which paramater(s) is(are) **required** to run the pipeline? 
- What is the default output directory's name? 
- What happens if you do not specify a profile (`-profile`)?

3.3 In the [nextflow-io *awesome* pipelines](https://github.com/nextflow-io/awesome-nextflow), look for the featured `BABS-aDNASeq` workflow:
- What tool is used for calling variants?
- What version of Nextflow is it advised to use?
- How do you download the `BABS-aDNASeq` pipeline locally?
````
````{tab} Solution 3
3.1. As of 20/10/2022: 71 pipelines are available, of which 40 are released, 25 are under development, and 6 are archived. 

3.2 [link](https://nf-co.re/atacseq)   
 - 1.2.2 (20/10/2022)  
 - 5 versions: current (1.2.2), 1.2.1, 1.2.0, 1.1.0, and 1.0.0.  
 - Only one required parameter: `--input` (Path to comma-separated file containing information about the samples in the experiment)  
 - `./results` (parameter `--outdir`)  
 - If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the PATH. More information is available [here](https://nf-co.re/atacseq/1.2.2/usage#main-arguments).   

3.3 [link](https://github.com/crickbabs/BABS-aDNASeq).  
 - `samtools mpileup`  
 - version 0.30.2 (Note that the current version is 22.10.0 (20/10/2022))  
 - `git clone https://github.com/crickbabs/BABS-aDNASeq`  (or `nextflow clone`)  

````
