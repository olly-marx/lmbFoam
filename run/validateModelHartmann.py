# Python script to vary the Ha number and run the model
# Then postprocess the results and save the data in the validation folder

import os
import sys

# Loop over the Ha number 1, 5, 20
Ha = [1, 5, 20]
# Loop over nYCells 10, 20, 40, 80, 160
nCells = [40]

# Take the solver name from the command line
solverNames = ["mhdFoam", "lmbFoam"]

# For each solver, run the model
for sol in solverNames:
    # For each Ha number, run the model
    for i in range(len(Ha)):
        # For each nCells, run the model
        for j in range(len(nCells)):
            dt = 2/(nCells[j]*1.2)

            #Echo the parameters
            print("Ha = " + str(Ha[i]))
            print("nCells = " + str(nCells[j]))
            print("dt = " + str(dt))
            print("solver = " + str(sol))

            os.system('bash ./testHartmannforHa_.sh ' + \
                    str(Ha[i]) + ' ' + str(nCells[j]) + ' ' + sol + ' ' + str(dt))

            # Postprocess the results
            os.system('pvpython ../../validation/processHartmannData.py ' + \
                    str(Ha[i]) + ' ' + str(nCells[j]) + ' ' + sol)
