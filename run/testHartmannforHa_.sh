#!/bin/sh

# Bash script used to run the hartmann case for a chosen value of Ha
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Get value of Ha from command line
Ha=$1

# Get the resolution from the command line
ny=$2

# Get the solver from the command line
solver=$3

# Get the stable time step from the command line
dt=$4

# Change the 0/Bext file line 19 to the chosen value of B_0
# If the solver is lmbFoam then edit Bext
# If the solver is mhdFoam then edit B
if [ $solver = "lmbFoam" ]
then
    dir=hartmannLMB
    sed -i "19s/.*/B0 $Ha;/" $dir/0/Bext
elif [ $solver = "mhdFoam" ]
then
    dir=hartmannMHD
    sed -i "19s/.*/B0 $Ha;/" $dir/0/B
fi

# Clear the case
pyFoamClearCase.py $dir

# Calculate the number of cells in x and y directions, based on the resolution
# ny = res, nx = 2.5*ny
nx=$(($ny*5/2))

# Change the blockMeshDict file to the chosen resolution, lines 31,32,33
sed -i "31s/.*/nx $nx;/" $dir/system/blockMeshDict
sed -i "32s/.*/ny $ny;/" $dir/system/blockMeshDict
sed -i "33s/.*/nz 1;/" $dir/system/blockMeshDict

# Change the controlDict file to the chosen time step, line 27
sed -i "27s/.*/deltaT $dt;/" $dir/system/controlDict

# Echo the simulation starting and the chosen values of Ha and total number of cells
nCells=$(($nx*$ny))
echo "Starting simulation for Ha = $Ha and $nx x $ny = $nCells cells"

cd $dir

# Run blockMesh
runApplication blockMesh

# Run the simulation
runApplication setFields
runApplication $solver

# And the .log files
mv log.$solver ./logs/log.Ha$Ha.N$ny.$solver
mv log.blockMesh ./logs/log.Ha$Ha.N$ny.blockMesh
mv log.setFields ./logs/log.Ha$Ha.N$ny.setFields

# paraFoam
