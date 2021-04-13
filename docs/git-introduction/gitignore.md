
# Gitignore
What if we have files that we do not want Git to track for us, like backup files or intermediate files created during data analysis? Remember that GitHub is not your next cloud storage infrastructure. Hence, your (big) data should not be uploaded on GitHub. Besides, there's a strict file size limit of 100MB so you won't even be able to do so. 

Regardless of the above, it is often useful if your data is in the same projectfolder. And you can't help the fact that Jupyter Notebooks makes intermediate checkpoints (.ipynb_checkpoints) in the same folder of the notebook. 

Git has a file, the `.gitignore` file, that has a list of expressions which define the files it doesn't need to track. This chapter will briefly discuss the `.gitignore` file without any exercises and going into detail. 

## Expressions
Imagine the following project folder structure:

```
 project-folder/
    |
    |- .ipynb_checkpoints
    |- .Rhistory
    |
    |- data/
    |   |- R1.fastq
    |   |- dataset.csv
    |
    ...
```

- **Ignore a file**:

The easiest would be to define the file or the path to the file. E.g. the fastq file can be ignored by adding `data/R1.fastq` to the `.gitignore` file. 

Similar to a file, a folder can also be ignored. The folders data/, .ipynb_checkpoints and .Rhistory can be ignored by adding the following lines:
```
data/
.ipynb_checkpoints
.Rhistory
``` 

- **`*`, `!` and `#`**:

The asterisk is often used in `.gitignore` files and represents a wildcard. E.g. `*.csv` will ignore any csv file in your folder. The asterisk can precede a file format in which case it will ignore all the files with that format (csv, fastq, sam, bam, xlsx, pdf, etc.) 

An exclamation mark is used for exceptions. The following lines of code will ignore all files in the data folder, except for the csv dataset:
```
data/
!data/dataset.csv
```

Documentation lines are preceded by a `#`. 

## Standard files

It's always good to think this through and manually add the files or folders that need to be ignored. However, it's also useful to know that there are standardized `.gitignore` files depending on the programming environment. They are all accessible in [this repository](https://github.com/github/gitignore) and contain `.gitignore` files for Python, R, Ruby, Java, Perl, C++, amongst many others. These files can also be added on the fly to a new repository by initializing the repository with one of these files (see figure below). 


```{image} ../img/git/gitignore.PNG
:align: center
```
