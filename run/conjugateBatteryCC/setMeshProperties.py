# Python script to find the system/meshProperties file and set the
# number of cells in the x, y, and z directions.
# Will take in a command line argument for the cell lateral size ratio
# and an optional argument to run the case in high resolution mode.

import os
import sys
import math

# Check for command line arguments
if len(sys.argv) < 2:
    print("Usage: python setMeshProperties.py <Lx/Ly ratio> [-HighRes (optional)]")
    sys.exit(1)
elif len(sys.argv) == 2:
    highRes = False
    LxLyRatio = float(sys.argv[1])
elif len(sys.argv) == 3:
    if sys.argv[2] == "-HighRes":
        highRes = True
        LxLyRatio = float(sys.argv[1])
    else:
        print("Usage: python setMeshProperties.py <Lx/Ly ratio> [-HighRes (optional)]")
        print("Error: Invalid argument: " + sys.argv[2])
        sys.exit(1)
else:
    print("Usage: python setMeshProperties.py <Lx/Ly ratio> [-HighRes (optional)]")
    sys.exit(1)

# Find the system/meshProperties file
meshPropertiesFile = "./system/meshProperties"

# Set the number of cells in the x, y, and z directions
N = 36

if highRes: Nx = 2*N; Ny = 2*N; Nz = 100
else: Nx = N; Ny = N; Nz = 50

NB = math.floor(N/8)

xyArea = 144*144
xyRatio = LxLyRatio

x1 = -0.5 * (xyArea / xyRatio)**0.5
x2 = 0.5 * (xyArea / xyRatio)**0.5
y1 = x1 * xyRatio
y2 = x2 * xyRatio
z1 = 0.0
z2 = 50.0

xLen1 = x2 - x1
yLen1 = y2 - y1
zLen1 = z2 - z1

deltax = xLen1/Nx
deltay = yLen1/Ny
deltaz = zLen1/Nz

print("delta x = " + str(deltax))
print("delta y = " + str(deltay))
print("delta z = " + str(deltaz))

x3 = -(Nx-2*NB)*deltax/2
x4 = (Nx-2*NB)*deltax/2
y3 = -(Ny-2*NB)*deltay/2
y4 = (Ny-2*NB)*deltay/2
z3 = 35.0
z4 = 50.0

xC1 = Nx-2*NB
xC2 = NB
yC1 = Ny-2*NB
yC2 = NB
zC1 = (z3-z1)/deltaz
zC2 = (z4-z3)/deltaz

x1b = x1 - 2.0
x2b = x2 + 2.0
y1b = y1 - 2.0
y2b = y2 + 2.0
z1b = z1 - 2.0
z2b = z2

x3b = x1
x4b = x2
y3b = y1
y4b = y2
z3b = z1
z4b = z2

xC1b = Nx
xC2b = 4
yC1b = Ny
yC2b = 4
zC1b = Nz
zC2b = 4

# Write the using sed commands
# The first line is the LxLyRatio, which needs to have 2 decimal places
os.system("sed -i 's/^LxLy .*/LxLy %.2f;/' %s" % (LxLyRatio, meshPropertiesFile))

# Now the rest
os.system("sed -i 's/^x1 .*/x1 %f;\/\/ \"$x1\"/' %s" % (x1, meshPropertiesFile)) 
os.system("sed -i 's/^x2 .*/x2 %f;\/\/ \"$x2\"/' %s" % (x2, meshPropertiesFile))
os.system("sed -i 's/^y1 .*/y1 %f;\/\/ \"$y1\"/' %s" % (y1, meshPropertiesFile))
os.system("sed -i 's/^y2 .*/y2 %f;\/\/ \"$y2\"/' %s" % (y2, meshPropertiesFile))
os.system("sed -i 's/^z1 .*/z1 %f;\/\/ \"$z1\"/' %s" % (z1, meshPropertiesFile))
os.system("sed -i 's/^z2 .*/z2 %f;\/\/ \"$z2\"/' %s" % (z2, meshPropertiesFile))
os.system("sed -i 's/^x3 .*/x3 %f;\/\/ \"$x3\"/' %s" % (x3, meshPropertiesFile))
os.system("sed -i 's/^x4 .*/x4 %f;\/\/ \"$x4\"/' %s" % (x4, meshPropertiesFile))
os.system("sed -i 's/^y3 .*/y3 %f;\/\/ \"$y3\"/' %s" % (y3, meshPropertiesFile))
os.system("sed -i 's/^y4 .*/y4 %f;\/\/ \"$y4\"/' %s" % (y4, meshPropertiesFile))
os.system("sed -i 's/^z3 .*/z3 %f;\/\/ \"$z3\"/' %s" % (z3, meshPropertiesFile))
os.system("sed -i 's/^z4 .*/z4 %f;\/\/ \"$z4\"/' %s" % (z4, meshPropertiesFile))
os.system("sed -i 's/^xC1 .*/xC1 %.0f;\/\/ \"$xC1\"/' %s" % (xC1, meshPropertiesFile))
os.system("sed -i 's/^xC2 .*/xC2 %.0f;\/\/ \"$xC2\"/' %s" % (xC2, meshPropertiesFile))
os.system("sed -i 's/^yC1 .*/yC1 %.0f;\/\/ \"$yC1\"/' %s" % (yC1, meshPropertiesFile))
os.system("sed -i 's/^yC2 .*/yC2 %.0f;\/\/ \"$yC2\"/' %s" % (yC2, meshPropertiesFile))
os.system("sed -i 's/^zC1 .*/zC1 %.0f;\/\/ \"$zC1\"/' %s" % (zC1, meshPropertiesFile))
os.system("sed -i 's/^zC2 .*/zC2 %.0f;\/\/ \"$zC2\"/' %s" % (zC2, meshPropertiesFile))
os.system("sed -i 's/^x1b .*/x1b %f;\/\/ \"$x1b\"/' %s" % (x1b, meshPropertiesFile))
os.system("sed -i 's/^x2b .*/x2b %f;\/\/ \"$x2b\"/' %s" % (x2b, meshPropertiesFile))
os.system("sed -i 's/^y1b .*/y1b %f;\/\/ \"$y1b\"/' %s" % (y1b, meshPropertiesFile))
os.system("sed -i 's/^y2b .*/y2b %f;\/\/ \"$y2b\"/' %s" % (y2b, meshPropertiesFile))
os.system("sed -i 's/^z1b .*/z1b %f;\/\/ \"$z1b\"/' %s" % (z1b, meshPropertiesFile))
os.system("sed -i 's/^z2b .*/z2b %f;\/\/ \"$z2b\"/' %s" % (z2b, meshPropertiesFile))
os.system("sed -i 's/^x3b .*/x3b %f;\/\/ \"$x3b\"/' %s" % (x3b, meshPropertiesFile))
os.system("sed -i 's/^x4b .*/x4b %f;\/\/ \"$x4b\"/' %s" % (x4b, meshPropertiesFile))
os.system("sed -i 's/^y3b .*/y3b %f;\/\/ \"$y3b\"/' %s" % (y3b, meshPropertiesFile))
os.system("sed -i 's/^y4b .*/y4b %f;\/\/ \"$y4b\"/' %s" % (y4b, meshPropertiesFile))
os.system("sed -i 's/^z3b .*/z3b %f;\/\/ \"$z3b\"/' %s" % (z3b, meshPropertiesFile))
os.system("sed -i 's/^z4b .*/z4b %f;\/\/ \"$z4b\"/' %s" % (z4b, meshPropertiesFile))
os.system("sed -i 's/^xC1b .*/xC1b %.0f;\/\/ \"$xC1b\"/' %s" % (xC1b, meshPropertiesFile))
os.system("sed -i 's/^xC2b .*/xC2b %.0f;\/\/ \"$xC2b\"/' %s" % (xC2b, meshPropertiesFile))
os.system("sed -i 's/^yC1b .*/yC1b %.0f;\/\/ \"$yC1b\"/' %s" % (yC1b, meshPropertiesFile))
os.system("sed -i 's/^yC2b .*/yC2b %.0f;\/\/ \"$yC2b\"/' %s" % (yC2b, meshPropertiesFile))
os.system("sed -i 's/^zC1b .*/zC1b %.0f;\/\/ \"$zC1b\"/' %s" % (zC1b, meshPropertiesFile))
os.system("sed -i 's/^zC2b .*/zC2b %.0f;\/\/ \"$zC2b\"/' %s" % (zC2b, meshPropertiesFile))

# Write the delta values to ./placeAttachmentPoints.py
placeAttachmentPointsFile = "./placeAttachmentPoints.py"

os.system("sed -i 's/^dx =  .*/dx =  %f/' %s" % (deltax, placeAttachmentPointsFile))
os.system("sed -i 's/^dy =  .*/dy =  %f/' %s" % (deltay, placeAttachmentPointsFile))
os.system("sed -i 's/^dz =  .*/dz =  %f/' %s" % (deltaz, placeAttachmentPointsFile))
