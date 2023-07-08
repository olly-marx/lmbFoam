### This script is used to run the parametric study for the liquid metal battery case
### It will first clean the case then subsequently run the case for each parameter combination
import os
import subprocess
import sys

# take in the parameter values from the terminal
# the first argument is the parameter value
# the second argument is the parameter name
value = float(sys.argv[1])
name = str(sys.argv[2])
boundaryConditions = str(sys.argv[3])
print("Parameter value: "+str(value))
print("Parameter name: "+str(name))
print("Boundary conditions: "+str(boundaryConditions))

### import the simple module from the paraview
from paraview.simple import *

# echo start of process in terminal 
print("START: Generating d2 data for parameter "+str(name)+" with value "+str(value)+" and boundary conditions "+str(boundaryConditions)+"...")

# create a new 'OpenFOAMReader'
simpleBatteryfoam = OpenFOAMReader(registrationName='simpleBattery.foam', FileName='./simpleBattery.foam')
simpleBatteryfoam.MeshRegions = ['internalMesh']
simpleBatteryfoam.CellArrays = ['alpha1']

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=simpleBatteryfoam)
contour1.ContourBy = ['POINTS', 'alpha1']
contour1.Isosurfaces = [0.5]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=contour1)
calculator1.Function = '0.034-coordsZ'

# create a new 'Descriptive Statistics'
descriptiveStatistics1 = DescriptiveStatistics(registrationName='DescriptiveStatistics1', Input=calculator1,
    ModelInput=None)
descriptiveStatistics1.VariablesofInterest = ['Result']

# set active source
SetActiveSource(descriptiveStatistics1)

# create save directory
os.system("mkdir -p ./postprocessing/d2Data"+str(name)+"/"+str(boundaryConditions))

# echo saving data in terminal
print("Saving d2 data...")

# save data
SaveData('./postprocessing/d2Data'+str(name)+'/d2.csv', proxy=descriptiveStatistics1, WriteTimeSteps=1,
    ChooseArraysToWrite=1,
    RowDataArrays=['Minimum'],
    FieldAssociation='Row Data',
    AddTime=1)

os.system("head -1 postprocessing/d2Data"+str(name)+"/d2_0.csv > postprocessing/d2All.csv")                  ## writing the header to the final file
os.system("tail -q -n +2  postprocessing/d2Data"+str(name)+"/d2_*.csv >> postprocessing/d2All.csv")             ## writing the content of all csv starting with second line into final file
os.system("rm postprocessing/d2Data"+str(name)+"/d2_*.csv")                                                     ## removing all csv files
os.system("mv postprocessing/d2All.csv postprocessing/d2Data"+str(name)+"/"+str(boundaryConditions)+"/d2Result"+str(name)+"_"+str(value)+".csv")           ## renaming the final file
    
# echo end of process in terminal
print("END: d2 data generated.")
    
