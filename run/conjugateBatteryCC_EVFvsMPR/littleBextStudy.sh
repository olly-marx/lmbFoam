#!/bin/bash

for Bext in 0.0005 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01
do
	# Change the Bext in the 0/Bext file
	sed -i "19s/internalField   uniform (0 0 -.*);/internalField   uniform (0 0 -$Bext);/" 0/Bext

	echo "Running with External B Field = $Bext"

	# Run the simulation
	./Allrun -quiet -parallel

	# Post process the simulation
	pvpython ../postprocessing/uRMS.py conjugateBatteryCC u_rms_Bext_$Bext

	# clean the case
	pyFoamClearCase.py .
	rm -rf processor* log*
done
