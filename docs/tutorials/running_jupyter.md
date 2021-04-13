
# Running Jupyter on Google Cloud

This tutorial will explain how to run a [Jupyter](https://jupyter.org/) Notebook or Jupyter Lab on a [Google Cloud](https://cloud.google.com/) instance. This set-up combines the intuitive and interactive assets of Notebooks along with markdown text documentation together with the strength of scaleable cloud computing provided by Google. 

Prerequisites: 
- Google Cloud instance up and running.

## 1. Installations

We'll have to start with some installations first:
- Google Cloud settings
- Installing Miniconda
- Installing Jupyter and related packages in a Conda environment
- Create the Jupyter configuration file

**Google Cloud settings** 

We will be running Jupyter in a local browser that is connected to the server. Hence, we need to make sure that this connection can be made. Therefore we need to make sure that Firewallrules are set properly. These steps are more elaborately discussed [here](https://towardsdatascience.com/running-jupyter-notebook-in-google-cloud-platform-in-15-min-61e16da34d52):
1. In Google Cloud, first find Compute Engine and click on VM-instance (at least if your instance is running there). Click on the instance from where you want to run Jupyter, click on 'Edit' and under Firewalls make sure that HTTP and HTTPS traffic are enabled. 
2. Go to VPC-network in the navigation menu and click on External IP-addresses and make a static IP-address for the server. You'll find your external IP-address now also under External IP-addresses between your VM-instances (e.g. 12.123.45.456)
3. Go to VPC-network in the navigation menu and click on Firewall(rules) and check that ports tcp:22 and tcp:3389 are enabled. These ports will enable SSH and RDP access on your instance.  


**Install Miniconda**  

Login to your Google Cloud instance and install the correct version of [Miniconda](https://docs.conda.io/en/latest/miniconda.html). Most probably you'll be interested in installing the 64-bit Linux installer with the latest Python version:

```none
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
# Press enter, type yes, check path is correct, yes to init
exec bash
```

After downloading Miniconda with `wget` you can verify that it was properly downloaded by verifying its hash with the following command `sha256sum <filename>` in which <filename> is the name of the downloaded file. It should be the same as the SHA256 hash given in the table [here](https://docs.conda.io/en/latest/miniconda.html).  

**Installing Jupyter within a conda environment**  
Now that we have (Mini)conda, we can easily make an environment with software packages that we might need for further data analysis. In this case we will make an environment with only Jupyter in it. 

```
conda activate
conda create -n running-jupyter jupyter
```

This means that first Conda is activated. Then, an environment is created using Conda with the name (`-n`) running-jupyter and with the package(s) jupyter. Basically you can extend this command with a list of packages that you want to install and thus use in the environment. 

To be verified useful packages.
- `nb_conda_kernels`
- `jupyterlab` - allows to run Jupyter Lab 
- `nodejs` 


**Jupyter configuration file**  
Finally, we need to check a couple of things in the Jupyter configuration file. 

You can set a password with the following command followed by entering a password:

```
jupyter notebook password
```

When logging in to the Notebook, Jupyter will ask you to confirm your login with this password. This will make a `jupyter_notebook_configuration.json` file with your encrypted password. 

Next, if it isn't there yet, you have to make the [configuration file](https://jupyter-notebook.readthedocs.io/en/stable/config_overview.html) which will make a `jupyter_notebook_configuration.py` file: 

```
jupyter notebook --generate-config
```

Check whether the following four items are in there [1](https://towardsdatascience.com/running-jupyter-notebook-in-google-cloud-platform-in-15-min-61e16da34d52):

```
c = get_config()
c.NotebookApp.ip = '*'
c.NotebookApp.open_browser = False
c.NotebookApp.port = <Port Number>
```

in which you replace `<Port Number>` by the port number which is enabled in the server and you will use later to connect to the Notebook (e.g. 5000 (or what else? why does it only work with htis number?))

## 2. Running Jupyter
You're now set to run a Jupyter Lab. In any directory and with the activated Conda environment, start a Jupyter Lab by:

``` 
jupyterlab --no-browser --port=<Port Number>
```

With `<Port Number>` e.g. 5000. 

Then open your browser (e.g. Google Chrome or Firefox) and navigate to: 

```
<external ip-address>:<port number>
```

With `<external ip-address>` in this case e.g. 12.123.45.456. And `<port number>` in this case e.g. 5000. 


Now, you might be interested to keep the notebook/lab running in the background[2](https://github.com/GenomicsCoreLeuven/vsc_ngs_workshop). For this part we can use either `screen` or `tmux`. Here we will elaborate on the latter. 

If you don't have Tmux installed, you can either follow their [installation instructions](https://github.com/tmux/tmux/wiki) or simply install it in your Conda environment:

```
conda install -c conda-forge tmux
```

Perform the following commands: 

```
# Open a tmux session
tmux
# Fire up a Jupyter notebook/lab 
jupyter lab --no-browser --port=<Port Number>
# Close the tmux session with 'Ctrl + b', followed by typing 'd'
# List your active tmux sessions with
tmux ls
```

To end session:

```
tmux kill-session -t <name-of-the-session>
```

In the example above, `<name-of-the-session>` will probably be 0. 

To reattach to the session which you detached from:

```
tmux attach -t <name-of-the-session>
```

## 3. References and further reading
- [1](https://towardsdatascience.com/running-jupyter-notebook-in-google-cloud-platform-in-15-min-61e16da34d52)
- [2](https://github.com/GenomicsCoreLeuven/vsc_ngs_workshop)
- [3](https://ljvmiranda921.github.io/notebook/2018/01/31/running-a-jupyter-notebook/)