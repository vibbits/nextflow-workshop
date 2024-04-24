```{toctree}
:hidden:
```

# Installations

Please read this page carefully **before** the start of the workshop.

There are two options for following this workshop: (1) do the installations yourself & be in control of everything, (2) use the setup that we have provided with the installations already done. In the former case, you will have to download [Nextflow](https://www.nextflow.io/docs/stable/getstarted.html) and [Apptainer](https://apptainer.org/). In the latter case, you can follow the instructions below.

## Provided infrastructure

We will be using the Gent section of the [Flemish Supercomputing Center](https://www.vscentrum.be/), you should have already recieved instructions for creating an account. Specifically, we will be using the [Interactive and Debug](https://docs.hpc.ugent.be/Linux/interactive_debug/) cluster. The cluster is already equipped with the latest version of Nextflow, and Apptainer.

To connect to the cluster, there are two options. For the first option, there is no local setup needed, we will use the Web Interface to access the Gent VSC.

### Option 1: Web Interface

This utilizes the OnDemand infrastructure at the VSC to launch a web-based version of VSCode for us. Using this, we don't need to make any connections to the clutser other than through the browser.

_If you normally use VSCode locally, this setup is completely seperate and won't have your usual extensions etc._

- Navigate to [https://login.hpc.ugent.be/](https://login.hpc.ugent.be/) and login with your credentials.
- Select "Interactive Apps" from the top bar -> "Code Server"
- Fill in the following settings:
  - Cluster: `donphan (interactive/debug)`
  - Time: 8 (hours)
  - Nodes: 1
  - Cores: 8
  - Select Path -> $VSC_DATA (on the left)
  - Click "Launch"!
- Wait for your job to start -> "Connect to VS Code"

This runs fully in your browser and will continue to run even when your laptop is off etc. Your job will automatically end after 8 hours. **Make sure to save your work.**

### Option 2: Local Installation

We will be using an SSH connection in VSCode which we can create by following these instructions:

- Download Visual Studio Code ([link](https://code.visualstudio.com/download))
- Add the following extensions for a seamless integration of Nextflow and the VM in VScode:
  - In VSCode, navigate to the 'Extensions' tab, search for the SSH remote package and install it:
  - 'Remote - SSH' (ms-vscode-remote.remote-ssh).
- Modify your local `.ssh/config` file to add the configuration for the cluster - If you already connect to the Gent VSC with this machine, you don't need to do this
  - `Ctrl-Shift-P` will bring up the "command palette"
  - Type `ssh config` and select the option to modify the configuration file (select the first file)
  - Add the following code to your config file:
    ```
    Host login-gent
        User vscXXXXX # Replace Xs with your VSC ID
        HostName login.hpc.ugent.be
        IdentityFile ~/.ssh/id_rsa # This should be replaced with the path to your private key ( windows users might look like this: C:\Users\KrisDavie\Documents\VSC\vsc_id_rsa)
    ```
- Start a terminal in VSCode (select Terminal and then New Terminal)
- Connect to the cluster with the following command: `ssh login-gent`
- Optional: Start `screen` or `tmux` and do the following in the new terminal - This will keep your session alive even when you disconnect from the cluster
- Load the modules for connecting to the interactive cluster: `module swap cluster/donphan`
- Start a job using qsub: `qsub -I -l walltime=08:00:00,nodes=1:ppn=8`
- Note the node you are connected to (e.g. `node4006.donphan.os`)
- Add the following code to your config file:
  ```
  Host node4006
      User vscXXXXX # Replace Xs with your VSC ID
      HostName node4006.donphan.os
      ProxyCommand ssh login-gent -W %h:%p
      # On windows you should use
      # ProxyCommand C:\Windows\System32\OpenSSH\ssh.exe login-gent -W %h:%p
  ```
- Finally you can open a this host in VSCode by typing `Ctrl-Shift-P` and selecting `Remote-SSH: Connect to Host...` and selecting the host you just added.
  - If you didn't run qsub in a screen or tmux session, you will need to use an entire new VSCode window to connect to the host, otherwise when VSCode refreshes, the original connection will be lost and the job will end.

### Option 3: Custom Installation

You are free to connect to the cluster however you want, but the above 2 methods are the only ones we will support in the session.

## Common Setup

- Install the Nextflow VSCcode Package - This will give you syntax highlighting and linting for Nextflow
- Open a new terminal within VSCode: Terminal -> New Terminal
- Create a new folder for the workshop
- Clone this repository into the folder: `git clone git@github.com:VIBbits/nextflow-workshop.git`
- For the Gent cluster usage, in any terminal where you want to run your excercises 
  1. Load the nextflow module: `module load Nextflow/23.10.0`
  2. Export the following envrionment variables - these are required so that your home folder is not filled when building and storing apptainer images
   - `export APPTAINER_CACHEDIR=${VSC_SCRATCH}/.apptainer_cache`
   - `export APPTAINER_TMPDIR=${VSC_SCRATCH}/.apptainer_tmp`
