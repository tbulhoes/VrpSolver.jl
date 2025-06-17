# The TOP application can be executed by invoking the docker directly:
docker run --rm -v /ABSOLUTE_PATH_TO_TOP_APP:/TOP bapdock /TOP/src/run.jl /TOP/data/p4.2.a.txt --cfg /TOP/config/TOP.cfg

# Interactive mode:
docker run -it --rm -v /ABSOLUTE_PATH_TO_TOP_APP:/TOP bapdock

# Help with command line arguments
docker run --rm -v /ABSOLUTE_PATH_TO_TOP_APP:/TOP bapdock /TOP/src/run.jl --help

# It is possible to run a batch of instances:
docker run --rm -v /ABSOLUTE_PATH_TO_TOP_APP:/TOP bapdock /TOP/src/run.jl -b /TOP/chao.batch

# If you are calling docker through a bash terminal (e.g. Linux, MacOS or Docker QuickStart Terminal), you can call the script named VRPSolver in the demo directory. For example:
./VRPSolver data/p4.2.a.txt --cfg config/TOP.cfg

# If you don't have permission to run VRPSolver script, call "chmod +x VRPSolver" before.
# This script must be called in the root directory of the application.

# Interactive mode:
./VRPSolver -it

# Help with command line arguments
./VRPSolver --help

# Running a batch of instances (see chao.batch before for adjustments):
./VRPSolver -b chao.batch

# Files with the extension .sh contain the call of VRPSolver for all instances individually.
