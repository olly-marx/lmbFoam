#!/bin/bash

Bext=$1

# Set the external magnetic field with sed in save/Bext
sed -i "s/internalField   (0 0 -.*/internalField   $Bext;/" save/Bext
sed -i "s/internalField   (0 0 -.*/internalField   $Bext;/" 0/Bext

for LxLy in 1.0 1.5 2.0 2.5 3.0 3.5 4.0
do
	# Set the cell dimensions with ./setMeshProperties.py
	python ../setMeshProperties.py $LxLy
	
	# Set the attachment point with ./placeAttachmentPoints.py
	python ../placeAttachmentPoints.py 0.0 0.0 -0.002

	# Run the simulation
	./Allrun -quiet -parallel

	# Post process the simulation
	pvpython ../postprocessing/ppData.py conjugateBatteryCC_Ratio_Bext RatioParaStudyBext-$LxLy-Bottom

	# clean the case
	pyFoamClearCase.py .
	rm -rf processor* log*
done
