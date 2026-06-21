# Legacy Installation (Docker)

Some advanced branching and cut separation functionalities are temporarily disabled in the current release of **VrpSolver.jl** (see [Disabled functions](disabled_functions.md)). If your workflow strictly requires those methods, you must use the **older Docker-based installation of VrpSolver**. This legacy environment contains the specific dependencies and core logic required for those routines to function.

The legacy Docker image (`bapdock.img`) and the demo applications built for the legacy VRPSolver are available here: [**Download the legacy files**](https://drive.google.com/drive/folders/1qMvajhSb2wJIq-lrcfXu-89qVtA4PI8O?usp=sharing).

The full step-by-step instructions are provided below.

## Installation instructions

This tutorial teaches you how to set up and run the legacy VRPSolver on Linux and on Windows/MacOS systems. While it assumes the use of CPLEX 12.10, it can be used in an equivalent way for installation with CPLEX 12.9.

```@raw html
<details>
<summary><strong>Linux</strong></summary>
```

First, follow the [instructions](https://docs.docker.com/engine/install/linux-postinstall/) to manage Docker as a non-root user (without `sudo`).

#### Installation of bapdock.img

Import and run the Docker image:

```bash
docker import -c "ENTRYPOINT [\"/julia\"]" bapdock.img bapdock-init
docker run -it -d --name bapdock-cont bapdock-init
```

Run this command to get the container-ID:

```bash
docker ps
```

Enter the subdirectory of your CPLEX installation which contains the shared libraries (e.g. `/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux`) and create a copy of `libcplex12100.so` with the name `libcplex12x.so`.

Copy your CPLEX folder to inside the container:

```bash
docker cp /ABSOLUTE_PATH_TO_CPLEX_ROOT/cplex bapdock-cont:/cplex
```

For example, if your CPLEX 12.10 is installed at `/opt/ibm/ILOG/CPLEX_Studio1210` this command should be:

```bash
docker cp /opt/ibm/ILOG/CPLEX_Studio1210/cplex bapdock-cont:/cplex
```

Commit the modified container as an image named `bapdock`:

```bash
docker commit container-ID bapdock
```

Stop and remove the `bapdock-init` container:

```bash
docker stop container-ID
docker rm container-ID
```

#### Running a Julia VRPSolver application through Docker

To run the application from your local machine command line, you must mount (`-v`) the user Julia VRPSolver application folder under any name, and the last arguments must be the main script `.jl` of the application, followed by the application arguments:

```bash
docker run --rm -v /ABSOLUTE_PATH_TO_VRPSOLVER_APP:/MyApp bapdock /MyApp/myapp.jl arg1 arg2 ... argn
```

The VRPSolver application folder is mounted as `/MyApp` in the Docker container filesystem. For example, for the *Capacitated Vehicle Routing Problem* (CVRP) demo (available among the [demos](https://drive.google.com/drive/folders/1qMvajhSb2wJIq-lrcfXu-89qVtA4PI8O?usp=sharing)), we can call:

```bash
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock /CVRP/src/run.jl /CVRP/data/A/A-n37-k6.vrp -m 6 -M 6 -u 950 -o /CVRP/sol/A-n37-k6.sol
```

which solves the instance `A-n37-k6.vrp` using exactly 6 vehicles (`-m` and `-M` define the minimum and maximum number of vehicles, respectively), an upper bound of 950, and writes the solution at `/CVRP/sol/A-n37-k6.sol`. Note that for the arguments after `bapdock` we need to consider the container filesystem (`/CVRP`) to access the application and I/O arguments. All changes into `/CVRP` and subdirectories will be reflected at `/ABSOLUTE_PATH_TO_CVRP_APP`.

When running from the local machine command line, Julia code will be compiled each time and thus building the model takes some time (like 30 seconds in a typical machine) even for small instances. The reported solution time, counted only after the code is compiled, is not affected.

To decrease recompilation time while creating/debugging the application, you can run the application in interactive mode:

```bash
docker run -it --rm -v /ABSOLUTE_PATH_TO_VRPSOLVER_APP:/MyApp bapdock
```

For example, for the CVRP demo, we can call:

```bash
docker run -it --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock
```

and to run the application inside the Julia environment, type (one at a time):

```julia
include("/CVRP/src/run.jl") # load CVRP demo
main(["/CVRP/data/A/A-n37-k6.vrp","-m","6","-M","6","-u","950","-o","/CVRP/sol/A-n37-k6.sol"]) # run CVRP demo
main(["/CVRP/data/A/A-n37-k5.vrp","-m","5","-M","5","-u","670"]) # running another instance without writing the solution
args = ["/CVRP/data/A/A-n37-k5.vrp","-m","5","-M","5","-u","671"] # defining arguments before
main(args) # running with predefined arguments
```

You can edit the application on your local machine and load it inside the Julia environment in Docker to test your changes, using `include`. In this case, recompilation time will be minimized.

To quit the Julia environment, type:

```julia
exit()
```

##### Easier way of running a Julia VRPSolver application

There is a bash script called `VRPSolver` available for each demo application to avoid calling Docker directly. For example, for the CVRP demo (see `README.txt` in the demo for more details), we can call:

```bash
./VRPSolver data/A/A-n37-k5.vrp -m 5 -M 5 -u 670
```

instead of:

```bash
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock /CVRP/src/run.jl /CVRP/data/A/A-n37-k5.vrp -m 5 -M 5 -u 670
```

```@raw html
</details>
```

```@raw html
<details>
<summary><strong>Windows or MacOS</strong></summary>
```

First, make sure Docker is in *running* status. If you are using Docker Toolbox, please use the *Docker QuickStart Terminal* to run Docker.

#### How to extract CPLEX

Create a folder named `cplexdocker`.

Get the CPLEX 12.10 installer file for Linux 64 from IBM ILOG (let's assume the file name is `cplex_studio1210.linux-x86-64.bin`) and move it to the folder `cplexdocker`.

Also, move `bapdock.img` to `cplexdocker`.

Create a text file named `Dockerfile_extract_cplex.txt` at the `cplexdocker` folder with this content:

```dockerfile
FROM ubuntu:24.04
RUN apt-get update ; apt-get install -y openjdk-8-jre
COPY cplex_studio1210.linux-x86-64.bin /
RUN (echo 2 ; echo ; echo 1 ; echo / ; echo y ; echo ; echo; echo; echo 2; echo ) | bash ./cplex_studio1210.linux-x86-64.bin
```

Run these commands at the directory `cplexdocker`:

```bash
docker build -t ubuntu_with_cplex -f Dockerfile_extract_cplex.txt .
docker create --name ubuntu_with_cplex ubuntu_with_cplex
```

Copy the installed CPLEX to the `cplexdocker` folder:

```bash
docker cp ubuntu_with_cplex:/cplex .
```

Enter the folder `cplex/bin/x86-64_linux/` and create a copy of `libcplex12100.so` with the name `libcplex12x.so`.

#### Installation of bapdock.img

Import and run the Docker image:

```bash
docker import -c "ENTRYPOINT [\"/julia\"]" bapdock.img bapdock-init
docker run -it -d --name bapdock-cont bapdock-init
```

Run this command to get the container-ID:

```bash
docker ps
```

Copy your CPLEX folder to inside the container:

```bash
docker cp cplex bapdock-cont:/cplex
```

Commit the modified container as an image named `bapdock`:

```bash
docker commit container-ID bapdock
```

Stop and remove the `bapdock-init` container:

```bash
docker stop container-ID
docker rm container-ID
```

#### Running a Julia VRPSolver application through Docker

To run the application from your local machine command line, you must mount (`-v`) the user Julia VRPSolver application folder under any name, and the last arguments must be the main script `.jl` of the application, followed by the application arguments:

```bash
docker run --rm -v /ABSOLUTE_PATH_TO_VRPSOLVER_APP:/MyApp bapdock /MyApp/myapp.jl arg1 arg2 ... argn
```

The VRPSolver application folder is mounted as `/MyApp` in the Docker container filesystem. For example, for the *Capacitated Vehicle Routing Problem* (CVRP) demo (available among the [demos](https://drive.google.com/drive/folders/1qMvajhSb2wJIq-lrcfXu-89qVtA4PI8O?usp=sharing)), we can call:

```bash
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock /CVRP/src/run.jl /CVRP/data/A/A-n37-k6.vrp -m 6 -M 6 -u 950 -o /CVRP/sol/A-n37-k6.sol
```

which solves the instance `A-n37-k6.vrp` using exactly 6 vehicles (`-m` and `-M` define the minimum and maximum number of vehicles, respectively), an upper bound of 950, and writes the solution at `/CVRP/sol/A-n37-k6.sol`. Note that for the arguments after `bapdock` we need to consider the container filesystem (`/CVRP`) to access the application and I/O arguments. All changes into `/CVRP` and subdirectories will be reflected at `/ABSOLUTE_PATH_TO_CVRP_APP`.

When running from the local machine command line, Julia code will be compiled each time and thus building the model takes some time (like 30 seconds in a typical machine) even for small instances. The reported solution time, counted only after the code is compiled, is not affected.

To decrease recompilation time while creating/debugging the application, you can run the application in interactive mode:

```bash
docker run -it --rm -v /ABSOLUTE_PATH_TO_VRPSOLVER_APP:/MyApp bapdock
```

For example, for the CVRP demo, we can call:

```bash
docker run -it --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock
```

and to run the application inside the Julia environment, type (one at a time):

```julia
include("/CVRP/src/run.jl") # load CVRP demo
main(["/CVRP/data/A/A-n37-k6.vrp","-m","6","-M","6","-u","950","-o","/CVRP/sol/A-n37-k6.sol"]) # run CVRP demo
main(["/CVRP/data/A/A-n37-k5.vrp","-m","5","-M","5","-u","670"]) # running another instance without writing the solution
args = ["/CVRP/data/A/A-n37-k5.vrp","-m","5","-M","5","-u","671"] # defining arguments before
main(args) # running with predefined arguments
```

You can edit the application on your local machine and load it inside the Julia environment in Docker to test your changes, using `include`. In this case, recompilation time will be minimized.

To quit the Julia environment, type:

```julia
exit()
```

##### Easier way of running a Julia VRPSolver application

If you are using the MacOS terminal or the *Docker QuickStart Terminal* (available for Docker Toolbox on Windows and MacOS), it is possible to use the bash script `VRPSolver` available for each demo application to avoid calling Docker directly. For example, for the CVRP demo (see `README.txt` in the demo for more details), we can call:

```bash
./VRPSolver data/A/A-n37-k5.vrp -m 5 -M 5 -u 670
```

instead of:

```bash
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock /CVRP/src/run.jl /CVRP/data/A/A-n37-k5.vrp -m 5 -M 5 -u 670
```

```@raw html
</details>
```

```@raw html
<details>
<summary><strong>Troubleshooting</strong></summary>
```

If you got the message `libcplex12x.so: cannot open shared object file: No such file or directory`, it means that the CPLEX directory was not successfully mounted. Please verify:

1. you created a copy of `libcplex12100.so` (or `libcplex1290.so`) named `libcplex12x.so`;
2. that the CPLEX version is 12.9 or 12.10;
3. that the CPLEX 64-bit Linux distribution is used, even if you have MacOS or Windows.

```@raw html
</details>
```
