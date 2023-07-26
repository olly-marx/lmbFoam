# Python script to set the attachement points for the battery
# in the OpenFOAM case
# The script takes three arguments for the position of the outer patch
# and the desired number of cells for the patch
# Usage: python placeAttachmentPoints.py <x> <y> <z> <N>

import os
import sys
import math

# Check for command line arguments
if len(sys.argv) < 4:
    print("Usage: python placeAttachmentPoints.py <x> <y> <z>")
    sys.exit()
elif len(sys.argv) == 4:
    x = float(sys.argv[1])
    y = float(sys.argv[2])
    z = float(sys.argv[3])
else:
    print("Usage: python placeAttachmentPoints.py <x> <y> <z>")
    sys.exit()

# The dx, dy, dz are constants and will be set by another script
dx =  5.000000
dy =  5.000000
dz =  1.000000

# Reset dx, dy, dz in based on N, and based on the value of z
#if z == -0.002 and y == 0.077:
#    dx=2.1
#    dy=2.1
#    dz=2.1
#elif z == -0.002:
#    dx=2.1
#    dy=2.1
#    dz=2.1
#elif z > -0.002:
#    dx=2.1
#    dy=2.1
#    dz=2.1

dx=3.1
dy=3.1
dz=3.1

# Scale the dx, dy, dz to m
dx = dx*0.001
dy = dy*0.001
dz = dz*0.001

x1 = x - dx
x2 = x + dx
y1 = y - dy
y2 = y + dy
z1 = z - dz
z2 = z + dz

# Write the delta values to system/splitPatchDict
splitPatchDictFile = "system/splitPatchDict"

os.system("sed -i 's/^x1  .*/x1  %f;/' %s" % (x1, splitPatchDictFile))
os.system("sed -i 's/^x2  .*/x2  %f;/' %s" % (x2, splitPatchDictFile))
os.system("sed -i 's/^y1  .*/y1  %f;/' %s" % (y1, splitPatchDictFile))
os.system("sed -i 's/^y2  .*/y2  %f;/' %s" % (y2, splitPatchDictFile))
os.system("sed -i 's/^z1  .*/z1  %f;/' %s" % (z1, splitPatchDictFile))
os.system("sed -i 's/^z2  .*/z2  %f;/' %s" % (z2, splitPatchDictFile))

# Now set the x, y, z values for the innerAttachmentPoint
x = 0.00
y = 0.00
z = 0.05

# Recalculate x1, x2, y1, y2, z1, z2
x1 = x - dx
x2 = x + dx
y1 = y - dy
y2 = y + dy
z1 = z - dz
z2 = z + dz

# Write the delta values to system/splitPatchDict2
splitPatchDictFile = "system/splitPatchDict2"

os.system("sed -i 's/^x1  .*/x1  %f;/' %s" % (x1, splitPatchDictFile))
os.system("sed -i 's/^x2  .*/x2  %f;/' %s" % (x2, splitPatchDictFile))
os.system("sed -i 's/^y1  .*/y1  %f;/' %s" % (y1, splitPatchDictFile))
os.system("sed -i 's/^y2  .*/y2  %f;/' %s" % (y2, splitPatchDictFile))
os.system("sed -i 's/^z1  .*/z1  %f;/' %s" % (z1, splitPatchDictFile))
os.system("sed -i 's/^z2  .*/z2  %f;/' %s" % (z2, splitPatchDictFile))
