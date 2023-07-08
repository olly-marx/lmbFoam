#!/bin/sh

# bash to loop over all the Bext runs in the data folder, restore them with the
# restoreRun.sh script and then run the processRunDataD1.py script to process the
# data 

for Bext in 0.0005 #0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0055 0.006 0.0065 0.007 0.0075 0.008 0.0085 0.009 0.0095 0.01
do
	# restore the case
	./restoreRun.sh data/Bext_$Bext

	# Run the python script to process the data
	python ./processRunDataD1.py $Bext Bext
	
	# now move the files into a new directory for safe keeping
	#rm -rf data/Bext_$Bext
	#mkdir -p data/Bext_$Bext
	#mv ./[0-9].* data/Bext_$Bext

	# clean the case
	#pyFoamClearCase.py .
done

