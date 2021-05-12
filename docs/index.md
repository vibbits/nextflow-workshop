---
hide-toc: true
---

# Containers & Workflow pipelines workshop

Welcome to the Containers & Workflow pipelines workshop! We are very happy to have you here. 


## General context
The first day (20 May 2021) is dedicated to Containers (Docker & Singularity) which are great tools for code portability and reproducibility of your analysis. You will learn how to use containers and how to build a container from scratch, share it with others and how to re-use and modify existing containers. After an extensive explanation on Docker containers, at the end of the first day, Singularity will be highlighted as well. 

On the second day (27 May 2021), you will learn how to use Nextflow for building scalable and reproducible bioinformatics pipelines and running them on a personal computer, cluster and cloud. Starting from the basic concepts we will build our own simple pipeline and add new features with every step, all in the new DSL2 language.  

## Objectives

- Containers
    - Learn the concept of and the difference between Docker & Singularity containers 
    - Pull and push Docker container to / from Dockerhub
    - Build and run Docker images and containers
    - Working with volumes (data inside your containers)
    - Understand Dockerfiles and layers; Docker cashing
    - Write a Docker recipe 


- Pipelines
    - Understand Nextflow's basic concepts: channels, processes, modules, workflows, etc. 
    - Write and run a Nextflow pipeline 
    - Write and modify config files for storing parameters related to computing hardware as well as pipeline dependent parameters

## Prerequisites 
- Being comfortable working with the CLI (command-line interface) in a Linux-based environment.

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
```

