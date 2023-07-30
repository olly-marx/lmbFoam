#!/bin/bash

for Voltage in 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 0.02 0.03 0.04 0.05 0.06 0.07
do
	# Change the Bext in the 0/Bext file
	sed -i "19s/cellVoltage   -.*);/cellVoltage   -$Voltage);/" 0/solid/solidVoltage

	echo "Running with cellVoltage = $Voltage"

	# Run the simulation
	./Allrun -quiet -parallel

	# Post process the simulation
	pvpython ../postprocessing/uRMS.py conjugateBatteryCC_Voltage dataFull-Voltage-$Voltage-Bottom

	# clean the case
	pyFoamClearCase.py .
	rm -rf processor* log*
done
