#!/bin/bash

for Bext in 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128
do
	# Change the Bext in the 0/Bext file
	sed -i "19s/internalField   uniform (0 0 -.*);/internalField   uniform (0 0 -$Bext);/" 0/Bext

	echo "Running with External B Field = $Bext"

	# Run the simulation
	./Allrun -quiet -parallel

	# Post process the simulation
	pvpython ../postprocessing/ppData.py conjugateBatteryCC_Bext Bext-$Bext-Bottom

	# clean the case
	rm -rf processor* log* [0-9]*
	cp -r save/ 0/
done
