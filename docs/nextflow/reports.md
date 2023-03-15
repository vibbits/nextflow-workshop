# Creating reports 
Nextflow has an embedded function for reporting a various information about the resources needed by each job and the timing. Just by adding a parameter on run-time, different kinds of reports can be created. 


1. **Workflow report**

    After running the nextflow pipeline script with the option `-with-report`, find the html report in the folder from where you launched the pipeline. 

    ```bash
    nextflow run exercises/05_reports/RNAseq.nf -with-report -profile docker
    ```
    
    This report describes the usage of resources and job durations and gives an indication of bottlenecks and possible optimizations in the pipeline. 

2. **DAG** 
    
    Use the option `-with-dag` to create a visualization of the workflow. By default and without any arguments, it will create a `.dot`-file that contains a description of the workflow, however to get a visualization we need to use an extra argument (e.g. `rnaseq.html`). This visualization is a nice overview of the workflow processes and how they are chained together and can be especially useful as a starting point to unravel more complex pipelines.

    ```bash
    nextflow run exercises/05_reports/RNAseq.nf -with-dag rnaseq.html -profile docker
    ```

    ```{note}
    As of Nextflow 22.04, the DAG can also be output in mermaid format, more information can be found [here](https://www.nextflow.io/docs/latest/tracing.html#dag-visualisation).
    ```

3. **Tower**

    Adding the parameter `-with-tower` enables the Seqera Tower service and will output the reports to a browser-based platform. More about Tower below.

## Tower
The Tower service, supported and developed by Seqera Labs, allows to monitor the workloads from a browser. Pipelines can be deployed on any local, cluster or cloud environment using the intuitive *launchpad* interface. Futhermore, it is also possible to manage teams and organizations, control project costs, and more. With ongoing improvements to the Tower platform, it is a very powerful platform worth checking out. 

To start using Tower, first create an account on [tower.nf](https://tower.nf). Then, we need to set the access token in our environment (on our VMs): 
```
export TOWER_ACCESS_TOKEN=<YOUR ACCESS TOKEN>
export NXF_VER=22.10.7 
```
Verify the Nextflow version (NXF_VER) with `nextflow -v`. The access token can be obtained from clicking on the top-right profile icon, select *Your tokens* and create *New token*. 

Tower is undergoing a lot of changes, hence we refer to this useful video. More information is also available at [seqera.io](https://seqera.io/).

<div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
    <iframe width="1280" height="720" src="https://www.youtube.com/embed/P7LUtBFzSww" title="Nextflow Tower" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
</div>

---


````{tab} Exercise 1
Run the `RNAseq.nf` pipeline again, this time also make the reports (both html-report and a visualization of the pipeline)
```` 
````{tab} Solution 1
The command that we need for this is the following.
```bash
nextflow run exercises/05_reports/RNAseq.nf -profile docker -with-report -with-dag rnaseq.html 
```
To view the report and the dag, you will need to download the files to your local machine.
```` 
