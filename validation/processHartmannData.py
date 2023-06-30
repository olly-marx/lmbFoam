# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *

# get active source.
hartmannfoam = GetActiveSource()

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=hartmannfoam)
plotOverLine1.Point1 = [10.0, -1.0, 0.05]
plotOverLine1.Point2 = [10.0, 1.0, 0.01]

# Properties modified on plotOverLine1Display_1
plotOverLine1Display_1.SeriesVisibility = ['B_X', 'U_X']

# get value of M from the command line
import sys
M = float(sys.argv[1])

# save data
SaveData('/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/validation/hartmann_lmbFoam_M' + str(M) + '.csv',
    proxy=plotOverLine1, ChooseArraysToWrite=1,
    PointDataArrays=['B', 'U'],
    AddTime=1)
