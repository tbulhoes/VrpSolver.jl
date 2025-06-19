# The GAP application can be executed by invoking the docker directly:
docker run --rm -v /ABSOLUTE_PATH_TO_GAP_APP:/GAP bapdock /GAP/src/run.jl /GAP/data/Classic/gapD-5-100.txt -u 6354 --cfg /GAP/config/GAP_Classic.cfg 

# toy instance (optimal is 449.0)
docker run --rm -v /ABSOLUTE_PATH_TO_GAP_APP:/GAP bapdock /GAP/src/run.jl /GAP/data/toy.txt 

# Interactive mode:
docker run -it --rm -v /ABSOLUTE_PATH_TO_GAP_APP:/GAP bapdock

# Help with command line arguments
docker run --rm -v /ABSOLUTE_PATH_TO_GAP_APP:/GAP bapdock /GAP/src/run.jl --help

# It is possible to run a batch of instances:
docker run --rm -v /ABSOLUTE_PATH_TO_GAP_APP:/GAP bapdock /GAP/src/run.jl -b /GAP/classic.batch

# If you are calling docker through a bash terminal (e.g. Linux, MacOS or Docker QuickStart Terminal), you can call the script named VRPSolver in the demo directory. For example:
./VRPSolver data/Classic/gapD-5-100.txt -u 6354 --cfg config/GAP_Classic.cfg

# If you don't have permission to run VRPSolver script, call "chmod +x VRPSolver" before.
# This script must be called in the root directory of the application.

# Interactive mode:
./VRPSolver -it

# Help with command line arguments
./VRPSolver --help

# Running a batch of instances (see classic.batch before for adjustments):
./VRPSolver -b classic.batch

# Files with the extension .sh contain the call of VRPSolver for all instances individually.
