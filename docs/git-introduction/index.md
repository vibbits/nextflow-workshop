# Git & GitHub

Git is an open-source tool that keeps track of changes to your files avoiding horrible situations depicted in the figure below. It works similar to [Google Docs'](https://support.google.com/drive/answer/2409045?co=GENIE.Platform%3DDesktop&hl=en) history feature in which Google automatically saves your document and the changes that happened on a particular moment in time. However, Git allows you to control and decide yourself when changes are worth saving, hence making it much more powerful and flexible. Each change is saved together with a message that enables you or your collaborators to keep an overview of the history of the project.  

```{toctree}
:hidden:

configurations
first_commit
history_status
branches
forks
gitignore
rstudio
```


```{image} ../img/git/version-control-meme.png
:align: center
```


Git is an open-source tool that manages your project (files) by keeping track of each version of these files throughout its history. It starts with a base version


Why should you version control? 
- **Keeping track of changes** to your files done by yourself or your collaborators. At any moment you can exploit the history of the project to see who wrote what on a particular day. It even allows you to go back to a specific version or undo specific edits. 
- **Synchronizes files between different people or infrastructures** (i.e. laptops, servers, ...), making it a powerful collaborating system. 
- **Testing new code**. Git can control multiple side versions of the same project in which you can make some changes and only when you or your collaborators are happy with hem, you can include them in the main version.


There is a major difference between Git and GitHub though. Git works on the command line of e.g. your computer, whereas GitHub is a service for connecting and uploading/downloading files much like saving files in the cloud. Alternatives for GitHub are Gitlab, Bitbucket, etc. In this course we will learn how Git works on the core of your computer which will give us proper understanding of its functionalities. Grasping these concepts is important if we want to use Git's version controlling in other apps (cfr. [GitHub and RStudio - insert link])

## Installations 
For this course we will explore version controlling in a mixture of [GitHub online](https://github.com/) & [Git](https://git-scm.com/) on the command line. The latter requires some basic understanding of the Linux command line. If you're not familiar with Linux command line, you can have a look at the materials [here](https://material.bits.vib.be/topics/linux/). After discussing Git's essential features, we'll introduce how you can make life easier with Git(Hub)'s integration in Rstudio. 

- Git can be installed for any OS (Windows, Mac or Linux) from [this link](https://git-scm.com/downloads).  
- Make an account on [GitHub](https://github.com/). 
- GitHub Desktop can be installed from [this link](https://desktop.github.com/). 

We will address further configurations in the next chapter. 

## Three conceptual areas
Before diving in, let's have a look at how Git works. It's important to understand the three conceptual areas that exist locally when using Git on your computer: the development area, the staging area and the repository containing the commits. 


```{image} ../img/git/conceptual_areas.png
:align: center
```


1. The **development area** is where your coding happens. Usually this is a folder with multiple files on your computer. Git will never change anything at this level, actually it won't really do anything. The only thing Git does is remembering that it needs to keep track of changes made in this folder or its files. Initializing Git on a project is only done once in the beginning. 
2. The **staging area** is an intermediate stage which assembles the files with changes you want to save. Imagine that we want to save a file, we first have to add it to the staging area before we can commit it.  
3. Files that are in the staging area are then committed to what we'll call the **commit repository** here. Committing is a synonym for saving the files in the Git terminology. The repository with commits contains a list of all the commits that we have done in a project. It's neatly structured in a history log which we can call at any point. Notice that all of this is still happening on our computer. 


Here's an example. Let's assume that we're starting a new project. Usually that  also means that you make a new folder on your computer where you will keep all the files related to the project. The first thing you have to do is to tell Git (or GitHub Desktop) that it has to keep track of this folder.In this step, we're initializing Git on this folder. Now, you just made your first file, however it's not automatically saved in Git. First, you'll have to add it to the staging area and afterwards you need to commit it to the repository. Voila, it's now safely stored in Git's repo.  
If we make a second file, the only thing we have to do now is to add it to the staging area and afterwards commit it. 

Notice that the repository is not yet visible on [github.com](https://github.com/). For this we would still need a fourth and last step, namely pushing the commits repository from your local machine to GitHub. By pushing your commits repository, you'll not only push the commits repository GitHub, but also push the files within the project to GitHub. After this last step, your project and all of the files are accessible in a GitHub repository.


During our adventure through git & GitHub we'll use some specific glossary. Confused on what the meaning of all these new words are? Check out the [GitHub glossary!](https://help.github.com/en/github/getting-started-with-github/github-glossary)