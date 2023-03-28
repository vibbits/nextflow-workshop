---
hide-toc: true
---

# Containers & Workflow pipelines workshop

Welcome to our Nextflow workshop! We are very happy to have you here. 


## General context
This is the third edition of this workshop, jointly organised by the VIB Bioinformatics Core and ELIXIR Belgium. 
- The first session (13 & 14 March 2023) is dedicated to Containers (Docker & Singularity) which are great tools for code portability and reproducibility of your analysis. You will learn how to use containers and how to build a container from scratch, share it with others and how to re-use and modify existing containers. 
- The second session (30 & 31 March 2023) is focused on Nextflow for building scalable and reproducible bioinformatics pipelines and running them on a personal computer, cluster and cloud. Starting from the basic concepts we will build our own simple pipeline and add new features with every step, all in the new DSL2 language. On the second day, we will utilise all the gathered knowledge to build a small-scale microbiomics pipeline. 

This website contains the course materials and outline for the second session. 

The presentation which goes alonside this material can be found [here](https://docs.google.com/presentation/d/1dl7yuVZTKeOKJwXuwTLb1NGWSZKKT0-THyllVtXMFsg/edit?usp=sharing).

## Practical information
Schedule day 1:

- 9:30 - 11:00 - session
- 11:00 - 11:15 - break
- 11:15 - 12:45 - session
- 12:45 - 13:45 - lunch
- 13:45 - 15:15 - session
- 15:15 - 15:30 - break
- 15:30 - 17:00 - session

*We aim to complete up to and including exercise 2.5 during this day*

Schedule day 2:

- 9:30 - 11:00 - session
- 11:00 - 11:15 - break
- 11:15 - 12:45 - session
- 12:45 - 13:45 - lunch
- 13:45 - 17:00 - project

## Objectives
The objectives of the Nextflow workshop are the following:
- Understand Nextflow's basic concepts & syntax: channels, processes, modules, workflows, etc. 
- Execute local and publicly available pipelines with different executors and environments 
- Write and run Nextflow pipelines  
- Write and modify config files for storing parameters related to computing hardware as well as pipeline dependent parameters

## Prerequisites 
Being comfortable working with the CLI (command-line interface) in a Linux-based environment.

## Requirements
The (technical) installation requirements are described in the [installations](https://nextflow-workshop.readthedocs.io/en/latest/installations.html) section. 

## Exercises
The exercises and solutions are available in [this GitHub repository](https://github.com/vibbits/nextflow-workshop).

```{toctree}
:hidden:

installations

```


```{toctree}
:hidden:
:caption: Courses

docker/index
nextflow/index
```



```{toctree}
:caption: Acknowledgements
:hidden:

ELIXIR Belgium <https://www.elixir-belgium.org/>
VIB Bioinformatics Core <https://www.bits.vib.be/>
BioLizard <https://www.lizard.bio>
```

