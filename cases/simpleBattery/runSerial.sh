#!/bin/sh

# bash for loop over ratios to run with ./runParamentricStudy.py script
# use args to give the script the ratio

for ratio in 1.0
do
	# Change the ratio in the blockMeshDict file
	sed -i "19s/.*atio.*/atio $ratio;/g" constant/polymesh/blockMeshDict

	# Set up blockMesh for this run
	source /usr/lib/openfoam/openfoam2212/etc/bashrc
	blockMesh
	source /home/ojm40/foam/foam-extend-5.0/etc/bashrc

	# Set up the initial conditions
	setFields

	# Run the program in series
	myLmbFoam

	# Run the python script to process the data
	python ./processRunDataD1.py $ratio
	
	# now move the files into a new directory for safe keeping
	#mkdir -p data/cellRatio_$ratio
	#mv ./*.* data/cellRatio_$ratio

	# clean the case
	#pyFoamClearCase.py .
done

# Run the python script to plot the data
python ./postprocessing/plotd1data.py
