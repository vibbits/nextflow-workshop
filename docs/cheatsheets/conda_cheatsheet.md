```{toctree}
:hidden:
```

# Conda cheatsheet

Activate conda base environment

```
$ conda activate
```

List all environments:

```
$ conda info --envs
```

Activate a specific environement (`environment-name`):

```
$ conda activate <environment-name>
```

Create a new environment:

```
$ conda create -n <environment-name> 
$ conda create -n <environment-name> <package>=<version> ... 
$ conda create -n <environment-name> -f <filename>.yml
```

In order to install a package from a channel in the current environment:

```
$ conda install -c <channel> <package>
```

Repository with packages: [Anaconda.org](Anaconda.org)


Export an environment (into a `yml`-file)
- In order to export current environment:

   ```
   $ conda env export > <filename>.yml
   ```

- Or, to export any other environment:

   ```
   $ conda env export -n <environment-name> > <filename>.yml
   ```

Remove a specific package from the current environment:

```
$ conda remove <package>
```

## Config file `.condarc` 
Show configuration file

```
$ conda config --show
```

Add channels

```
$ conda config --add channels conda-forge 
$ conda config --add channels bioconda
```

List channels

```
$ conda config --get channels
```

Add a channel with lowest priority

```
$ conda config --append channels newchannel
```

Add a channel with highest priority

```
$ conda config --prepend channels newchannel
```
