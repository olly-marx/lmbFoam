#!/bin/sh

# bash for loop over ratios 1.0, 1.5, 2.0 to run with ./runParamentricStudy.py script
# use args to give the script the ratio
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# take in args and use those to build the specific case to be run
# then run it in parallel and plot the output
# args are: Voltage, Bext, ratio

# Take args from the command line
Param=$1
Setting=$2
Voltage=$3

# Set defaults for Bext, Ratio and Voltage
sed -i "19s/internalField   uniform (0 0 -.*);/internalField   uniform (0 0 -0.0005);/" 0/Bext

# Check which parameter is being set in the run change the appropriate file and run the case
if [ $Param = "Bext" ];
then
	sed -i "19s/internalField   uniform (0 0 -.*);/internalField   uniform (0 0 -$Setting);/" 0/Bext
	Ratio=1.0
	Bext=$Setting
elif [ $Param = "Ratio" ];
then
	sed -i "19s/atio.*/atio $Setting;/g" constant/polyMesh/blockMeshDict
	Bext=0.0005
	Ratio=$Setting
fi

if [ $Voltage = "Sidewalls" ];
then
	# Active the voltage boundary condictions in the 0/Voltage file
	echo "Setting sidewalls"
	# Set the side wall to fixed value
	sed -i "25s/.*/\/\/\ttype\t\tzeroGradient;/" 0/Voltage
	sed -i "26s/.*/\ttype\t\tfixedValue;/" 0/Voltage
	sed -i "27s/.*/\tvalue\t\tuniform 0;/" 0/Voltage
	# Set the bottom wall to zero gradient
	sed -i "31s/.*/\ttype\t\tzeroGradient;/" 0/Voltage
	sed -i "32s/.*/\t\/\/type\t\fixedValue;/" 0/Voltage
	sed -i "33s/.*/\t\/\/value\t\tuniform 0;/" 0/Voltage
elif [ $Voltage = "Bottom" ];
then
	# Active the voltage boundary condictions in the 0/Voltage file
	echo "Setting bottom"
	# Set the side wall to zero gradient
	sed -i "25s/.*/\ttype\t\tzeroGradient;/" 0/Voltage
	sed -i "26s/.*/\t\/\/type\t\tfixedValue;/" 0/Voltage
	sed -i "27s/.*/\t\/\/value\t\tuniform 0;/" 0/Voltage
	# Set the bottom wall to fixed value
	sed -i "31s/.*/\t\/\/type\t\tzeroGradient;/" 0/Voltage
	sed -i "32s/.*/\ttype\t\tfixedValue;/" 0/Voltage
	sed -i "33s/.*/\tvalue\t\tuniform 0;/" 0/Voltage
fi

echo "Running with $Param = $Setting and Voltage = $Voltage, (Defaults) Bext = $Bext, Ratio = $Ratio, Current collector = $Voltage"

# Set up blockMesh for this run
source /usr/lib/openfoam/openfoam2212/etc/bashrc
blockMesh
source /home/ojm40/foam/foam-extend-5.0/etc/bashrc

# Set up the initial conditions
runApplication setFields

# contruct parallel mesh
runApplication decomposePar

# Run the program in parallel
mpirun -np 8 myLmbFoam -parallel

# Reconstruct the parallel data
runApplication reconstructPar

# Run the python script to process the data
python ./processRunDatad2.py $Setting $Param $Voltage

# count the number of time steps present and print to screen
echo "Number of time steps completed:"
find [0-9]* -maxdepth 0 | wc -l

# now move the files into a new directory for safe keeping
dir="$Param$Voltage$Setting"
echo "Moving files to $dir"

rm -rf data/$dir
mkdir -p data/$dir
mv ./[0-9]* data/$dir

mv ./log.* ./data/$dir/

# copy other case files to the data directory
cp -r ./system ./data/$dir/
cp -r ./constant ./data/$dir/
cp -r ./save ./data/$dir/

# now clear the case for the next run
pyFoamClearCase.py .
