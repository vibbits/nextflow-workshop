```{toctree}
:hidden:
```

# Linux cheatsheet
 
Directory size:

```
du -sh directory/
du -a -h --max-depth=1 | sort -hr
```

Server login:

```
ssh -i ssh/ssh-key user@11.111.11.111 
```

with `-i` the identity, and `ssh-key` private ssh key. 

Copy file from server to local 

```
scp -r -i ssh/ssh-key user@11.111.11.111:/home/directory/server/ /home/directory/local/
```

