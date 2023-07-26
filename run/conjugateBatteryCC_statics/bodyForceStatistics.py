# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
conjugateBatteryCCfoam = OpenFOAMReader(registrationName='conjugateBatteryCC.foam', FileName='/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/conjugateBatteryCC/conjugateBatteryCC.foam')
conjugateBatteryCCfoam.MeshRegions = ['internalMesh']
conjugateBatteryCCfoam.CellArrays = ['bodyForce']

UpdatePipeline(time=0.0012, proxy=conjugateBatteryCCfoam)

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=conjugateBatteryCCfoam)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.ResultArrayName = 'magBodyForce'
calculator1.Function = 'mag(bodyForce)'

UpdatePipeline(time=0.0012, proxy=calculator1)

# create a new 'Descriptive Statistics'
descriptiveStatistics1 = DescriptiveStatistics(registrationName='DescriptiveStatistics1', Input=calculator1,
    ModelInput=None)
descriptiveStatistics1.VariablesofInterest = ['magBodyForce']

UpdatePipeline(time=0.0012, proxy=descriptiveStatistics1)

UpdatePipeline(time=0.0012, proxy=descriptiveStatistics1)

# set active source
SetActiveSource(descriptiveStatistics1)

# save data
SaveData('/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/conjugateBatteryCC/results/bodyForceStats.csv', proxy=descriptiveStatistics1, ChooseArraysToWrite=1,
    RowDataArrays=['Maximum', 'Mean', 'Minimum', 'Variable'],
    Precision=8,
    FieldAssociation='Row Data',
    AddTime=1)
