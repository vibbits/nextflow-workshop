# Project

For day two, each of you will build a small-scale metagenomics pipeline from scratch to test out what you’ve learned so far.

## Concept

In metagenomics, environmental samples are taken and examined for the micro-organisms that can be found in them. In its most basic form, this is done by extracting DNA from each sample and selectively amplifying and sequencing the 16S rRNA-encoding gene, which is unique to bacteria and archaea. By examining the variety in sequences we get from sequencing only this gene, we can start to tell which specific micro-organisms are present within a sample.

In this project, we’d like to combine publicly available tools and some basic R scripts to:

1. Check the quality of the reads
2. Trim primers from them and filter low-quality reads
3. Find the unique 16S sequence variants & get a grasp of the diversity of the samples.

## Data

For this project, we will use data from *[Vandeputte et al. (2017)](https://www.nature.com/articles/nature24460)*: a well-known study from the VIB centre for microbiology that was published in Nature. This study took faecal samples from 135 participants to examine their gut microbiota.

To keep computation times to a minimum, we will work with a few subsets from this data. You can pull in this data by running the `get_project_data.sh` script that has been provided for you.

```bash
bash get_project_data.sh
```

## Our metagenomics pipeline

### Step 0: Preparation

We’d like to construct a pipeline that executes quality control, trimming, filtering and the finding of unique sequence in an automated fashion, parallelising processes wherever possible. We’d also like to run the whole thing in docker containers so we don’t have to worry about dependencies.

````{tab} Objective 1
Set up a `main.nf` script in which you will build your pipeline that uses nextflow’s DSL2 and which reads in the forward and reverse reads for each of the five samples in the data1-directory into a channel.
```` 
---

All the docker containers we will need are already publicly available, so don’t worry about having to write Dockerfiles yourself 🙂

### Step 1: Quality Control

After pulling in and setting up the data, we’re first interested in examining the quality of the sequencing data we got. 

````{tab} Objective 2 
Write a process which executes FastQC over the raw samples.
````
---

As we’re not really looking forward to inspecting each FastQC report individually, we should pool these in a single report using MultiQC.

````{tab} Objective 3
Write a second process that executes MultiQC on the directory of FastQC outputs.
```` 
---

If this all works, you should be able to take a look at the outputted `.html` report, in which you should see stats for 10 sets of reads (forward and reverse for each of the 5 samples). 

### Step 2: Trimming and filtering

Looking at the MultiQC report, our reads don’t look that fantastic at this point, so we should probably do something about that. 

We can use a publicly available tool called Cutadapt: 

- to trim off the primers
- to trim and filter low-quality reads
- to remove very short reads and reads containing unknown bases (*i.e.,* ‘N’)

Cutadapt however requires us to specify the forward and reverse primers, as well as their reverse [complements](http://www.reverse-complement.com/ambiguity.html). The forward and reverse primers we can find in the paper: `GTGCCAGCMGCCGCGGTAA` and `GGACTACHVHHHTWTCTAAT`, but the reverse complements we’ll have to figure out ourselves.

````{tab} Objective 4 
Put your nextflow operator knowledge to the test and write some nextflow code that calculates the reverse complements of these primers.
```` 
---

Once we have all of the primer sequences, we can run cutadapt on all the reads. In bash, the code for this would look something like this:

```bash
cutadapt -a ^FW_PRIMER...REVERSECOMP_RV_PRIMER \\
				 -A ^RV_PRIMER...REVERSECOMP_FW_PRIMER \\
				 --quality-cutoff 28 \\
				 --max-n 0 \\
				 --minimum-length 30 \\
				 --output SAMPLE_R1_TRIMMED.FASTQ --paired-output SAMPLE_R2_TRIMMED.FASTQ \\
				 SAMPLE_R1.FASTQ SAMPLE_R2.FASTQ
```

````{tab} Objective 5 
Write a process that executes Cutadapt to filter and trim the reads.
````
---

### Step 3: Re-evaluate

As hopefully Cutadapt has done job, we’d now like to take another look at the quality report of the preprocessed reads to see if this has improved the stats.

````{tab} Objective 6
Write a workflow in your `[main.nf](http://main.nf)` file which runs FastQC and MultiQC on the raw reads, filters and trims these reads using Cutadapt, and then reruns FastQC and MultiQC on the preprocessed reads.
````
````{tab} Hint
Combine the FastQC and MultiQC processes into a named workflow.
```` 
--- 
### Step 4: Find unique sequences and plot

To closely examine amplicon sequencing data and to extract from these the unique 16S sequence variants, there is an incredibly useful package in R called DADA2. You have been provided with a small R script which uses this package to count the abundance of each unique sequence in each sample. Moreover, this script also outputs a dendrogram of a hierarchical clustering of the samples based on this abundance table, giving us an idea of the beta-diversity of our samples. This script takes the preprocessed forward & reverse reads (in no specific order) as input arguments on the command line.

````{tab} Objective 7 
Write and incorporate a process that executes this Rscript and outputs the `counts_matrix.csv` and `dendrogram.png` files.
````
````{tab} Hint
The container that you use should have the R-package ‘DADA2’ installed.
````
---
You now have successfully written your own microbiomics pipeline!

### Step 5: Rinse and repeat

To see if all this effort in automatisation was really worth it, you should run your pipeline on another dataset to see if it works as well.

````{tab} Objective 8
Run your pipeline on the sequencing data in the `data2` directory.
```` 
---

### If you have time left

There are a few things left that you can implement in your pipeline so others can more easily work with it as well.

- Use directives to output data from different processes to separate directories.
- Make a cool header that displays every time you run your pipeline using the `[log.info](http://log.info)` command.
- Add an onComplete printout to your pipeline that tells the user where they can find the output files.
- Speed up the slow processes in your pipeline by allocating more cpus and memory to them.
- Have nextflow create a report when you run the pipeline to see some cool stats.
- Implement this pipeline in [tower.nf](http://tower.nf) and run it from the web-application