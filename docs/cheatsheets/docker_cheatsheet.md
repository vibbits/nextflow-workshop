```{toctree}
:hidden:
```

# Docker cheatsheet
 
Version of Docker:

```
$ docker --version
```

List all Docker containers

```
$ docker container ls
$ docker ps -a 
```

List all Docker images

```
$ docker image ls
```

Pull a Docker image

```
$ docker pull <image-name>
```

Repositories: [Quay.io](https://quay.io) & [Dockerhub](https://hub.docker.com/). The naming convention for Docker containers is: `OWNER/CONTAINER:TAG`. 


## Running parameters

Create and run container from Docker image

```
$ docker run <image-name>
```

- `-it`: interactively (use `exit` to finish using the container)
- `--rm`: remove container after running
`-v`: add volumes, i.e. using data within a container. Example: `-v ~/path/to/data/:/workdir` or for multiple volumes: `-v ~/path/to/data1/:/workdir/x -v ~/path/to/data2/:/workdir/y`

- `-w`: manually set the working directory (useful when using volumes). Example: `-w="/home"` 

- `--name`: label the container with a name 

- `-u $(id -u):$(id -g)`: output files created by the container carry the user- and group-levels (e.g. `user:group` instead of `root:root` 
- Entry level with e.g.: `/bin/bash -c "<command>"`


## Troubleshooting

1. Activating Docker 

    ```
    Connecting Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running? 
    ```
    Make sure that Docker is working. 

2. Docker root 
    ```
    Got permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock: Get http://%2Fvar%2Frun%2Fdocker.sock/v1.40/containers/json: dial unix /var/run/docker.sock: connect: permission denied
    ```
    Solve this issue with [this nice explanation](https://www.digitalocean.com/community/questions/how-to-fix-docker-got-permission-denied-while-trying-to-connect-to-the-docker-daemon-socket). 