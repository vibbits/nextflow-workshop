```{toctree}
:hidden:
```

# Installations 

Please read this page carefully before the start of the workshop. 

There are two options for following this workshop: (1) do the installations yourself & be in control of everything, (2) use the VMs that we have created with the installations already done. In the former case, you will have to download [Nextflow](https://www.nextflow.io/docs/edge/getstarted.html), [Docker](https://docs.docker.com/get-docker/) and [Singularity](https://sylabs.io/singularity/). In the latter case, you can follow the instructions below. 

## Provided infrastructure
We have created a classroom lab in Azure Lab for each participant. This means that each participant has access to a dedicated Azure Labs Virtual Machine, this is an Ubuntu 18.04 virtual machine containing all the software needed for the workshop, and which can be accessed via RDP or SSH. 
- RDP represents Remote Desktop Protocol and allows you to access the VM as a remote computer. 
- SSH enjoys our preference as it will be more stable and faster. It represents a Secure Shell Protocol and will be a command-line type of interaction with the VM.   
Each classroom lab can be accessed during the workshop hours and for an extra 10 hours. This extra time allows you to experiment with the materials in your spare time. 

To obtain the RDP and SSH keys, you can follow [this link](https://docs.microsoft.com/en-us/azure/lab-services/how-to-use-classroom-lab). Briefly, navigate to the registration URL that you received in your mailbox, make an account, find the VM in your Azure Labs, start the Lab and copy the RDP/SSH key. 

There are several ways of accessing the VM via SSH:
- On Linux and Mac via the Terminal
- On Windows you can install the Windows terminal or an application like MobaXterm.
**However**, we propose the integration with Visual Studio Code for all OS (Windows, Linux & Mac)! 
- Download Visual Studio Code ([link](https://code.visualstudio.com/download))
- In VSCode, navigate to the 'Extensions' tab, search for the following packages and install them:
    - 'Remote - SSH' (ms-vscode-remote.remote-ssh). 
    - Optionally: 'Nextflow' (nextflow.nextflow), 'Docker' (ms-azuretools.vscode-docker). 
- Navigate to the new 'Remote Explorer' icon, find the dropdown menu 'WSL Targets' and switch to 'SSH Targets'. Click on '+', in the pop-up navigation bar paste the SSH key & follow the instructions. 


## Communication platform
We will use Zoom for this training. Participants should have received the link to the Zoom meeting. 

