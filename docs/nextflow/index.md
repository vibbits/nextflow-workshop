

# Nextflow 
This tutorial aims to get you familiarized with Nextflow. After this course you should be able to understand workflow pipelines that are written in Nextflow and write simple pipelines yourself! Here's an overview of the materials that we will cover:

- General introduction to Nextflow 
- Basic concepts: processes, channels and operators
- Creating our first Nextflow script(s)
- Managing configurations: parameters, portability, execution
- Creating reports

```{toctree}
:hidden:

basic_concepts
first_pipeline
configs
reports
```

After a brief introduction and example in DSL1, we will immediately change to the newer version of Nextflow with DSL2. It's important to have an idea of the differences between the two versions as most of the pipelines today are written in DSL1, but being transformed or enriched with new pipelines written in DSL2. 

<!--
On the day of writing (November 2020), DSL2 has been introduced for quite a while and DSL1 is supposed to be fading out. Chances are that support for DSL1 will be gone within a year or so.)
--> 

## Prerequisites
This course requires familiarity with the command-line. You should feel confident interacting with the command-line and launch commands and tools from there. Additionally, Docker containers will be used during the course. 

## Installations
The following software should be installed on your machine:
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Docker](https://docs.docker.com/engine/install/)
- Optionally, [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

## References and further reading
Here are some great tips of where you can learn more and where to get inspiration for writing your own pipelines: 
- Nextflow's official documentation ([link](https://www.nextflow.io/docs/latest/index.html))
- Reach out to the community on Gitter ([link](https://gitter.im/nextflow-io/nextflow))
- Curated collection of patterns ([link](https://github.com/nextflow-io/patterns))
- Workshop focused on DSL2 developed by CRG Bioinformatics Core ([link](https://github.com/biocorecrg/ELIXIR_containers_nextflow))
- Tutorial exercises (DSL1) developed by Seqera ([link](https://github.com/seqeralabs/nextflow-tutorial))
- Curated ready-to-use analysis pipelines by NF-core ([link](https://nf-co.re/))
- Model example pipeline on Variant Calling Analysis with NGS RNA-Seq data developed by CRG ([link](https://github.com/CRG-CNAG/CalliNGS-NF))






