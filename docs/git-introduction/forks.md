
# Forking
In this chapter we will discuss a strategy for collaborating on a project. 

Imagine that you're starting a project with some colleagues and you want to version control the project. If it were to be a document where each of you needs to write part of it, you could simply start a Google Doc. For coding purposes the situation is a bit more complex. There might be a base version of the code already to which you need to add separate parts, however you always need to test whether your part is working together with the rest of the code. 

For this purpose, GitHub encourages the Fork & Pull workflow. Basically one **forks** a central repository, making it a personal forked repository. This repository can constantly be up to date with the central repository. Changes that are done by yourself are pushed to the forked personal repository. If the work is complete, the changes in the forked repository can be merged back into the central repository by creating a pull request. This workflow leaves the central repository untouched untill the moment you want to incorporate changes. 

Forking a remote repository is not enough. After you've forked a repository, it will appear as a new repository in your GitHub account. The next step would be to **clone** the repository locally so you can work on the project from your computer. It's always a good idea to make changes in a **new branch** and keep the *master* branch clean. Hence, after cloning the repository, you could make a new branch. Staging, committing and pushing your changes remains the same and they will appear in your new forked repository. 

```{image} ../img/git/fork_pull.png
:align: center
```



Two important terms in this fork & pull workflow are:
- `upstream`: generally refers to the original repository that you have forked
- `origin`: is your fork: your own repo on GitHub  

As mentioned in section 4.4, the "origin" is used to refer to the GitHub original repository's URL. This also lasts here. The remote `origin` refers to your fork on GitHub, not the original repository it was forked from. 

To summarize the above, the Fork & Pull workflow consists of the following steps:
1. Fork
2. Clone
3. Branch
4. Stage-commit-push
5. Pull request

## Fork
Let's first start with exploring a bit on GitHub. GitHub is like the Facebook of programmers. You can see someone's account, what that person has been working on, find new projects (relatable to a Facebook page), etc. Exploring new repositories is possible by clicking on the 'Explore' button in the navigation bar. Searching a specific repository or searching for an account, on the other hand, is possible by simply typing it in the search bar in the navigation bar. 


```{image} ../img/git/nav-bar.PNG
:align: center
```


Search for the VIB Bioinformatics Core account 'vibbits' and find the repository 'fork-repository'. This repository was made especially for learning the principle of forking. Do this by clicking the fork button in the upper right corner.



```{image} ../img/git/fork-button.PNG
:align: center
```

The repository has been successfully forked if you see something similar to the figure below. The icon represents a fork, followed by your GitHub account name and the name of the repository. Below it tells us that the upstream repository is the vibbits/forked-repository. 

```{image} ../img/git/forked-repository.PNG
:align: center
```

## Changes
In normal circumstances, this would be the point where you clone this repository locally, make a branch and do some edits in that branch. The flow here remains the same: stage-commit-push. For this exercise we will only edit the `participants.txt` file on GitHub as we've explored the local usage before. 

Add your name, accountname or initials and the date to the `participants.txt` file. 

Your repository now looks like this. Notice the indicator saying that this branch is 1 commit ahead of the upstream repository. 

```{image} ../img/git/edited-forked-repository.PNG
:align: center
```

This is of course only the case as there were no changes in the upstream repository in the meantime. In normal circumstances the upstream repository might have changed. The indicator would then note that there are new commits in the upstream (`1 commit behind vibbits:master`), while the branch/repo itself is one commit ahead.  

```{image} ../img/git/forked-repository-ahead.PNG
:align: center
```


This does not (really) affect the pull request.  

## Pull request
The two repositories have diverged during the previous steps. Now its time to create a pull request between these repositories. Find the pull request button right underneath the green clone or download button and click it. 

A new screen pops up that looks very similar to the one seen in Chapter 5 (Branching & merging). 


```{image} ../img/git/forked-pull-request.PNG
:align: center
```


GitHub tells us:
- It compared the master branch of the forked repository (in my case *tmuylder/fork-repository*) with the upstream (base) repository *vibbits/fork-repository*. 
- It's able to merge these two branches without any conflicting errors
- It summarizes the changes that have been done in the branch that will be merged into the upstream.  

If all seems good, you can create the pull request. In the case that there are any conflicting errors, they'll need to be solved first. Afterwards you only need to add a message that accompanies the pull request. 

A brief overview of the pull request is given in the following screen which either allows you to merge the pull request into the upstream repository yourself or which requests the maintainer of the upstream repository to review and merge the pull request. In the latter case, the maintainer will thereafter receive a notification showing the pull request. An overview of all pending pull requests where you are involved in, are consultable on the [pull requests](https://github.com/pulls) tab of the navigation bar.   


```{image} ../img/git/forked-repository-final-pull-request.PNG
:align: center
```

## Overview
This is the easiest collaboration you'll probably do in a lifetime. To briefly summarize, the steps that we took were: *fork > edit > pull request (> merge)*. As mentioned before this is only possible if the upstream repository didn't change (too much). If this were to be the case, there might be one additional step in which you have to solve conflicts in the pull request. 

If your changes were a bit more complex and needed to be performed on your local computer, the steps would extent to: *fork > clone (> branch) > edit-stage-commit-push > pull request (> merge)*. What if the upstream repository changed while you were working on your local repository? In this case a pull request should be done in which the receiving branch is your forked repository. Hence, the order of the branches as depicted in the figure above would be swapped.    



---

> ### {% icon hands_on %} Exercise 
>
> Merge upstream changes in your forked repository. This approach is useful if your working on a project that is prone to lots of changes and you need to keep up to date. 
> Note: This exercise is only possible to be performed if the repository `vibbits/fork-repository` has changed after you forked it.  
> 
>    > <details markdown="1">
>    > <summary>{% icon solution %} Solution
>    > </summary>
>    > You need to merge any upstream changes into your version, and you can do this with a pull request on GitHub too. This time though you will need to switch the bases of the  comparison around, because the changes will be coming from the upstream version to yours. First find the following notification in your repository and click on pull request:  
>    > <center><img src="../../images/Exercise-fork-1.PNG" /></center>
>    > In my case, the order is not how it's supposed to be and the message reads: "There isn't anything to compare. vibbits:master is up to date with all commits from tmuylder:master.". Click on *switching the base* in order to insert the changes from the upstream in your forked repository.  
>    > 
>    > A message similar to the following will allow to create a pull request and subsequently merge the changes into your forked repository. 
>    > 
>    > 
>    > <center><img src="../../images/Exercise-fork-2.PNG" /></center>
>    > 
>    > 
>    > </details>
> 
{: .hands_on}

---
 