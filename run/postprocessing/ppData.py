# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# import numpy
import numpy as np

# Take .foam file name as input
import sys
foamFile = sys.argv[1]

# get active source.
thisfoam = OpenFOAMReader(registrationName=foamFile+'.foam', \
        FileName='./'+foamFile+'.foam')
# Properties modified on thisfoam
thisfoam.MeshRegions = ['internalMesh']
thisfoam.CellArrays = ['U', 'alpha1', 'bodyForce']

# Contour the interface
contourI = Contour(registrationName='Contour1', Input=thisfoam)
contourI.PointMergeMethod = 'Uniform Binning'
contourI.ContourBy = ['POINTS', 'alpha1']
contourI.Isosurfaces = [0.5]

# Calculate height of interface
calculatorHi = Calculator(registrationName='CalculatorHI', Input=contourI)
calculatorHi.ResultArrayName = 'height_i'
calculatorHi.Function = 'coordsZ'

# Threshold the electrolyte
thresholdE = Threshold(registrationName='ThresholdE', Input=thisfoam)
thresholdE.Scalars = ['POINTS', 'alpha1']
thresholdE.UpperThreshold = 0.5

# Calculate height of electrolyte
calculatorHE = Calculator(registrationName='CalculatorHE', Input=thresholdE)
calculatorHE.ResultArrayName = 'height'
calculatorHE.Function = 'coordsZ'

# Threshold the part of the electrolye above the interface
thresholdErel = Threshold(registrationName='ThresholdErel', Input=calculatorHE)
thresholdErel.Scalars = ['POINTS', 'height']
thresholdErel.LowerThreshold = 0.026
thresholdErel.UpperThreshold = 0.036

# Calculate u rms in the electrolyte
calculatorUrmsInE = Calculator(registrationName='CalculatoruRMSinE',\
        Input=thresholdErel)
calculatorUrmsInE.ResultArrayName = 'magU2'
calculatorUrmsInE.Function = 'mag(U)^2'

# Plot over time all fluid
plotDataOverTime1 = PlotDataOverTime(registrationName='PlotDataOverTime1',\
        Input=thisfoam)

# Plot over time all fluid
plotDataOverTime2 = PlotDataOverTime(registrationName='PlotDataOverTime2',\
        Input=calculatorUrmsInE)

# Plot over time all fluid
plotDataOverTime3 = PlotDataOverTime(registrationName='PlotDataOverTime3',\
        Input=calculatorHi)

# set active source
SetActiveSource(plotDataOverTime1)

# get save name from input
saveName = sys.argv[2]

# set active source
SetActiveSource(plotDataOverTime1)

# save data
SaveData('../postprocessing/UvsTime.csv', proxy=plotDataOverTime1, ChooseArraysToWrite=1,
    RowDataArrays=['Time', 'avg(U (Magnitude))', 'avg(bodyForce (Magnitude))'],
    Precision=8,
    FieldAssociation='Row Data',
    AddTime=1)

# set active source
SetActiveSource(plotDataOverTime3)

# save data
SaveData('../postprocessing/URMSElecvsTime.csv', proxy=plotDataOverTime2, ChooseArraysToWrite=1,
    RowDataArrays=['Time','avg(magU2)'],
    Precision=8,
    FieldAssociation='Row Data',
    AddTime=1)

# set active source
SetActiveSource(plotDataOverTime2)

# save data
SaveData('../postprocessing/HvsTime.csv', proxy=plotDataOverTime3, ChooseArraysToWrite=1,
    RowDataArrays=['Time','max(height_i)'],
    Precision=8,
    FieldAssociation='Row Data',
    AddTime=1)

# Now take the previous data and combine it into one file, with the first column being time
# and the second column being the u_ave, and the third column being the
# bodyForce_ave, and the fourth column being the u_rms_elec, and the fifth column being the height_i
# get data from file ./UvsTime.csv ignoring the first row
data1 = \
np.genfromtxt('../postprocessing/UvsTime.csv',\
        delimiter=',')
# get data from file ./URMSElecvsTime.csv ignoring the first row
data2 = \
np.genfromtxt('../postprocessing/URMSElecvsTime.csv',\
        delimiter=',')
# get data from file ./HvsTime.csv ignoring the first row
data3 = \
np.genfromtxt('../postprocessing/HvsTime.csv',\
        delimiter=',')

# Now combine the data
data = np.column_stack((data1[:,0], data1[:,1], data1[:,2], data2[:,1], data3[:,1]))

# Now save the data
np.savetxt('../postprocessing/'\
        +saveName+'.csv', data, delimiter=',', header='Time, U_ave, bodyForce_ave, u_rms_elec, height_i')

