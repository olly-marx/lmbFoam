#!/bin/sh

# Bash script used to run the hartmann case for a chosen value of Ha
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Get value of Ha from command line
Ha=$1

# Get the resolution from the command line
ny=$2

# Get the solver from the command line
solver=$3

# Clear the case
pyFoamClearCase.py .

# And the .log files
rm log.*

# Change the 0/Bext file line 19 to the chosen value of B_0
# If the solver is lmbFoam then edit Bext
# If the solver is mhdFoam then edit B
if [ $solver = "lmbFoam" ]
then
    sed -i "19s/.*/B_0 $Ha;/" 0/Bext
    sed -i "19s/.*/B_0 $Ha;/" save/Bext
elif [ $solver = "mhdFoam" ]
then
    sed -i "19s/.*/B_0 $Ha;/" 0/B
    sed -i "19s/.*/B_0 $Ha;/" save/B
fi

# Calculate the number of cells in x and y directions, based on the resolution
# ny = res, nx = 2.5*ny
nx=$(($ny*5/2))

# Change the blockMeshDict file to the chosen resolution, lines 31,32,33
sed -i "31s/.*/nx $nx;/" system/blockMeshDict
sed -i "32s/.*/ny $ny;/" system/blockMeshDict
sed -i "33s/.*/nz 1;/" system/blockMeshDict

# Echo the simulation starting and the chosen values of Ha and total number of cells
nCells=$(($nx*$ny))
echo "Starting simulation for Ha = $Ha and $nx x $ny = $nCells cells"

# Run blockMesh
runApplication blockMesh

# Run the simulation
runApplication setFields
runApplication $solver
