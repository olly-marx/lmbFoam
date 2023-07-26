# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
conjugateBatteryfoam = OpenFOAMReader(registrationName='conjugateBattery.foam', FileName='/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/conjugateBattery/conjugateBattery.foam')
conjugateBatteryfoam.MeshRegions = ['internalMesh']
conjugateBatteryfoam.CellArrays = ['U', 'alpha1']

# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=conjugateBatteryfoam)
threshold1.Scalars = ['POINTS', 'alpha1']
threshold1.UpperThreshold = 0.5

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=threshold1)
calculator1.ResultArrayName = 'height'
calculator1.Function = 'coordsZ'

# create a new 'Threshold'
threshold2 = Threshold(registrationName='Threshold2', Input=calculator1)
threshold2.Scalars = ['POINTS', 'height']
threshold2.LowerThreshold = 0.026
threshold2.UpperThreshold = 0.034

# create a new 'Calculator'
calculator2 = Calculator(registrationName='Calculator2', Input=threshold2)
calculator2.ResultArrayName = 'u_rms'
calculator2.Function = 'mag(U)'

# create a new 'Plot Data Over Time'
plotDataOverTime1 = PlotDataOverTime(registrationName='PlotDataOverTime1', Input=calculator2)

# save data
SaveData('/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/conjugateBattery/results/uMagAve_vs_t.csv', proxy=plotDataOverTime1, ChooseArraysToWrite=1,
    RowDataArrays=['Time', 'avg(u_rms)'],
    Precision=8,
    FieldAssociation='Row Data',
    AddTime=1)