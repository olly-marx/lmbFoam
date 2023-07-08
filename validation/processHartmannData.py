# trace generated using paraview version 5.11.1
#import paraview

#### import the simple module from the paraview
from paraview.simple import *

# Get the solver from the command line
solver = sys.argv[3]
if(solver=="lmbFoam")
    directory = "hartmannLMB"
elif(solver=="mhdFoam")
    directory = "hartmannMHD"

# get active source.
# create a new 'OpenFOAMReader'
hartmannfoam = OpenFOAMReader(registrationName=directory+'.foam', \
FileName='/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/hartmann/'+directory+'.foam'
hartmannfoam.MeshRegions = ['internalMesh']
hartmannfoam.CellArrays = ['Bext', 'U']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on animationScene1
animationScene1.AnimationTime = 1

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=hartmannfoam)
plotOverLine1.Point1 = [10.0, -1.0, 0.05]
plotOverLine1.Point2 = [10.0, 1.0, 0.05]

# Get the Ha number from the command line
Ha = sys.argv[1]

# Get the NCells from the command line
NCells = sys.argv[2]

# Set the file name
fileName = '/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/validation/' + \
        'hartmann_' + solver + '_Ha' + Ha + '_NCells' + NCells + '.csv'

# save data
SaveData(fileName, proxy=plotOverLine1, ChooseArraysToWrite=1,
    PointDataArrays=['U'],
    AddTime=1)
