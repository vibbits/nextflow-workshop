```{toctree}
:hidden:
```

# Installations 

Please read this page carefully **before** the start of the workshop. 

There are two options for following this workshop: (1) do the installations yourself & be in control of everything, (2) use the VMs that we have created with the installations already done. In the former case, you will have to download [Nextflow](https://www.nextflow.io/docs/edge/getstarted.html), [Docker](https://docs.docker.com/get-docker/) and [Singularity](https://sylabs.io/singularity/). In the latter case, you can follow the instructions below. 

## Provided infrastructure
We have created a VM with Azure Labs for each participant. This means that each participant has access to a dedicated Azure Labs VM. This is an Ubuntu 18.04 virtual machine containing all the software needed for the workshop. Each VM can be accessed during the workshop hours and for an extra 10 hours. This extra time allows you to experiment with the materials in your spare time. 

There are several ways for connecting to the VMs, however we advise to connect with the VM through a SSH connection in VSCode following these instructions:
- Download Visual Studio Code ([link](https://code.visualstudio.com/download))
- Add the following extensions for a seamless integration of Nextflow and the VM in VScode:
    - In VSCode, navigate to the 'Extensions' tab, search for the following packages and install them:
    - 'Remote - SSH' (ms-vscode-remote.remote-ssh). 
    - 'Nextflow' (nextflow.nextflow) 
    - 'Docker' (ms-azuretools.vscode-docker). 
- Create the connection to the VM:
    - In VScode, navigate to the new 'Remote Explorer' icon on the left side bar, find the dropdown menu 'WSL Targets' and switch to 'SSH Targets'. 
    - Click on '+', in the pop-up navigation bar and paste the SSH key (see below) & follow the instructions. 
- To obtain the SSH keys:
    - Navigate to the registration URL of Microsoft Azure Lab that you received in your mailbox and select or make a Microsoft account
    - Find the VM *containers-workflow-2022* in your Azure Labs, start the Lab and choose a password. After some waiting, you can copy the SSH key.  
    - An elaborate description is available [here](https://docs.microsoft.com/en-us/azure/lab-services/how-to-use-classroom-lab). 

In a VScode Terminal (select Terminal and then New Terminal), test the infrastructure with the following command `nextflow -h`. The server will be empty, however we will populate it with the course materials during the workshop. 
