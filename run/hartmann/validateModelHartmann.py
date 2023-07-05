# Python script to vary the Ha number and run the model
# Then postprocess the results and save the data in the validation folder

import os
import sys

# Loop over the Ha number 1, 5, 20
Ha = [1, 5]
# Loop over nYCells 10, 20, 40, 80, 160
nCells = [10, 20, 40, 80, 160]

# Take the solver name from the command line
solverName = str(sys.argv[1])

# For each Ha number, run the model
for i in range(len(Ha)):
    # For each nCells, run the model
    for j in range(len(nCells)):
        os.system('bash ./testHartmannforHa_.sh ' + \
                str(Ha[i]) + ' ' + str(nCells[j]) + ' ' + solverName)

        # Postprocess the results
        os.system('pvpython ../../validation/processHartmannData.py ' + \
                str(Ha[i]) + ' ' + str(nCells[j]) + ' ' + solverName)
