# The CVRP application can be executed by invoking the docker directly:
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock /CVRP/src/run.jl /CVRP/data/A/A-n37-k6.vrp -m 6 -M 6 -u 950 -o /CVRP/sol/A-n37-k6.sol

# toy instance (optimal 265)
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock /CVRP/src/run.jl /CVRP/data/toy.vrp

# Interactive mode:
docker run -it --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock

# Help with command line arguments
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock /CVRP/src/run.jl --help

# It is possible to run a batch of instances:
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP bapdock /CVRP/src/run.jl -b /CVRP/A.batch

# The application directory (/ABSOLUTE_PATH_TO_CVRP_APP) was mounted with -v as /CVRP inside the container. Also, it is possible to mount a different directory to read/write solutions:
docker run --rm -v /ABSOLUTE_PATH_TO_CVRP_APP:/CVRP -v /ABSOLUTE_PATH_TO_OUTPUT:/OUT bapdock /CVRP/src/run.jl /CVRP/data/A/A-n37-k6.vrp -m 6 -M 6 -u 950 -o /OUT/A-n37-k6.sol

# If you are calling docker through a bash terminal (e.g. Linux, MacOS or Docker QuickStart Terminal), you can call the script named VRPSolver in the demo directory. For example:
./VRPSolver data/A/A-n37-k6.vrp -m 6 -M 6 -u 950 -o sol/A-n37-k6.sol

# If you don't have permission to run VRPSolver script, call "chmod +x VRPSolver" before.
# This script must be called in the root directory of the application.

# Interactive mode:
./VRPSolver -it

# Help with command line arguments
./VRPSolver --help

# Running a batch of instances (see A.batch before for adjustments):
./VRPSolver -b A.batch

# Files with the extension .sh contain the call of VRPSolver for all instances individually.
