# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# load state
LoadState('/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/conjugateBatteryCC/currentStreamsC.pvsm')

# USe command line argument to get the conductivity ratio
import sys
ratio = float(sys.argv[1])

# find source
plotOverLine2 = FindSource('PlotOverLine2')

# set active source
SetActiveSource(plotOverLine2)

# Properties modified on plotOverLine2
plotOverLine2.SamplingPattern = 'Sample At Segment Centers'

# save data
fileName = '/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/conjugateBatteryCC/results/solidCurrentAlongBase_' + str(ratio) + '.csv'
SaveData(fileName, proxy=plotOverLine2, ChooseArraysToWrite=1,
    PointDataArrays=['solidJ0'],
    Precision=8,
    AddTime=1)
