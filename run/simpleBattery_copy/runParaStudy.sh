#!/bin/sh

# bash for loop over ratios 1.0, 1.5, 2.0 to run with ./runParamentricStudy.py script
# use args to give the script the ratio
. $WM_PROJECT_DIR/bin/tools/RunFunctions

for Voltage in "Sidewalls" "Bottom"
do
	for Bext in 0.0005 0.005 0.05 0.5 1.0
	do
		# Change the Bext in the 0/Bext file
		sed -i "19s/internalField   uniform (0 0 -.*);/internalField   uniform (0 0 -$Bext);/" 0/Bext

		echo "Running with External B Field = $Bext and Boundary Conditions on $Voltage"

		# Set up blockMesh for this run
		source /usr/lib/openfoam/openfoam2212/etc/bashrc
		runApplication blockMesh
		source /home/ojm40/foam/foam-extend-5.0/etc/bashrc

		# Set up the initial conditions
		runApplication setFields

		# contruct parallel mesh
		runApplication decomposePar

		# Run the program in parallel
		runApplication mpirun -np 8 lmbFoam -parallel

		# Reconstruct the parallel data
		runApplication reconstructPar

		# Run the python script to process the data
		python ./processRunDatad2.py $Bext Bext $Voltage

		# count the number of time steps present and print to screen
		echo "Number of time steps completed:"
		find [0-9]* -maxdepth 0 | wc -l

		# now move the files into a new directory for safe keeping
		dir="Bext$Bext$Voltage"
		echo "Moving files to $dir"

		# Store the data using storeRun.sh
		bash ./storeRun.sh $dir

		# clean the case
		pyFoamClearCase.py .
	done

	# Run both the python scripts to plot the data
	python ./postprocessing/plotd2data.py Bext $Voltage
	python ./postprocessing/plotd2minimum.py Bext $Voltage

	# Reset the Bext in the 0/Bext file
	sed -i "19s/internalField   uniform (0 0 -.*);/internalField   uniform (0 0 -0.0005);/" 0/Bext

	#for ratio in 1.0
	#do
	#	# Change the ratio in the blockMeshDict file
	#	sed -i "19s/atio.*/atio $ratio;/g" constant/polyMesh/blockMeshDict

	#	echo "Running with ratio = $ratio and Boundary Conditions on $Voltage"

	#	# Set up blockMesh for this run
	#	source /usr/lib/openfoam/openfoam2212/etc/bashrc
	#	runApplication blockMesh
	#	source /home/ojm40/foam/foam-extend-5.0/etc/bashrc

	#	# Set up the initial conditions
	#	runApplication setFields

	#	# contruct parallel mesh
	#	runApplication decomposePar

	#	# Run the program in parallel
	#	runApplication mpirun -np 8 lmbFoam -parallel

	#	# Reconstruct the parallel data
	#	runApplication reconstructPar

	#	# Run the python script to process the data
	#	python ./processRunDatad2.py $ratio Ratio $Voltage

	#	# count the number of time steps present and print to screen
	#	echo "Number of time steps completed:"
	#	find [0-9]* -maxdepth 0 | wc -l

	#	
	#	# now move the files into a new directory for safe keeping
	#	dir="Ratio$ratio$Voltage"
	#	echo "Moving files to $dir"

	#	# Store the data using storeRun.sh
	#	bash ./storeRun.sh $dir

	#	# clean the case
	#	pyFoamClearCase.py .
	#done

	## Reset the ratio in the blockMeshDict file
	#sed -i "19s/atio.*/atio 1.0;/g" constant/polyMesh/blockMeshDict

	## Run the python script to plot the data
	#python ./postprocessing/plotd2data.py Ratio $Voltage
	#python ./postprocessing/plotd2minimum.py Ratio $Voltage

	# Active the voltage boundary condictions in the 0/Voltage file
	# Set the side wall to zero gradient
	sed -i "25s/.*/\ttype\t\tzeroGradient;/" 0/Voltage
	sed -i "26s/.*/\t\/\/type\t\tfixedValue;/" 0/Voltage
	sed -i "27s/.*/\t\/\/value\t\tuniform 0;/" 0/Voltage
	# Set the bottom wall to fixed value
	sed -i "31s/.*/\t\/\/type\t\tzeroGradient;/" 0/Voltage
	sed -i "32s/.*/\ttype\t\tfixedValue;/" 0/Voltage
	sed -i "33s/.*/\tvalue\t\tuniform 0;/" 0/Voltage

done

# Active the voltage boundary condictions in the 0/Voltage file
# Set the side wall to fixed value
sed -i "25s/.*/\/\/\ttype\t\tzeroGradient;/" 0/Voltage
sed -i "26s/.*/\ttype\t\tfixedValue;/" 0/Voltage
sed -i "27s/.*/\tvalue\t\tuniform 0;/" 0/Voltage
# Set the bottom wall to zero gradient
sed -i "31s/.*/\ttype\t\tzeroGradient;/" 0/Voltage
sed -i "32s/.*/\t\/\/type\t\fixedValue;/" 0/Voltage
sed -i "33s/.*/\t\/\/value\t\tuniform 0;/" 0/Voltage
