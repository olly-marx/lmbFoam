# Varying the position of the attachment point along a line, by running the
# placeAttachmentPoint.py script with different values of the x, y, z coordinates
# We will run "pyFoamClearCase.py ." followed by the script with a new set of
# coordinates, then run the case using ./Allrun, calculating and storing the
# maximum value of the body force for each run. We will then plot the maximum
# body force along this line in space.

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Define the line in space which goes from (0,0,-0.002) to (0,0.077,-0.002)
# and then from (0,0.077,-0.002) to (0,0.077,0.05)

# Define the number of points along the line
numPoints = 22

# Define an array to store the x, y, z coordinates of the points along the line
line1 = np.zeros((numPoints, 3))
line2 = np.zeros((numPoints, 3))

# Keep an array of the length travelled along the line
line1Length = np.zeros(numPoints)
line2Length = np.zeros(numPoints)

# Define 15 points along the first part of the line, all with x=0, z=-0.002,
# with y varying by 0.005 each time from 0 to 0.075
for i in range(13):
    line1[i,0] = 0.0
    line1[i,1] = 0.006*i
    line1[i,2] = -0.002
    line1Length[i] = 0.006*i

    line2[i,0] = 0.077
    line2[i,1] = 0.006*i
    line2[i,2] = -0.002
    line2Length[i] = 0.006*i

line1[13,0] = 0.0
line1[13,1] = 0.077
line1[13,2] = -0.002
line1Length[13] = 0.077

line2[13,0] = 0.077
line2[13,1] = 0.077
line2[13,2] = -0.002
line2Length[13] = 0.077

# Define 10 points along the second part of the line, all with x=0, y=0.077,
# with z varying by 0.005 each time from -0.002 to 0.048
for i in range(8):
    line1[i+14,0] = 0.0
    line1[i+14,1] = 0.077
    line1[i+14,2] = -0.002 + 0.006*(i+1)
    line1Length[i+14] = 0.077 + 0.006*(i+1)

    line2[i+14,0] = 0.077
    line2[i+14,1] = 0.077
    line2[i+14,2] = -0.002 + 0.006*(i+1)
    line2Length[i+14] = 0.077 + 0.006*(i+1)

## Now looping over the different points along the line, run the placeAttachmentPoint.py
## script, then run the case, then calculate and store the max, min and mean body force
## for each case, an array of 26 rows and 3 columns
#bodyForceStats = np.zeros((numPoints, 3))
#
#for i in range(numPoints):
#    # Run the placeAttachmentPoint.py script
#    print("Running placeAttachmentPoint.py with coordinates: "+str(line1[i,0])+ \
#            " "+str(line1[i,1])+" "+str(line1[i,2]))
#    os.system("python placeAttachmentPoints.py "+str(line1[i,0])+" "+ \
#            str(line1[i,1])+" "+str(line1[i,2]))
#    # Clear the case
#    os.system("pyFoamClearCase.py .")
#    # Run the case
#    os.system("./Allrun -special -serial")
#    # Calculate and store the maximum body force
#    os.system("pvpython bodyForceStatistics.py")
#    # Sanity check, if i=13, run paraFoam
#    if i == 13:
#        os.system("paraFoam")
#    elif i == 20:
#        os.system("paraFoam")
#    # Read in the data
#    data = np.genfromtxt("./results/bodyForceStats.csv", delimiter=",", \
#            skip_header=1)
#    # Store the maximum, minimum and mean body force, columns 1, 3 and 2
#    # Row 0
#    bodyForceStats[i,0] = data[1]
#    bodyForceStats[i,1] = data[3]
#    bodyForceStats[i,2] = data[2]
#    
## Save the data to a csv file
#np.savetxt("bodyForceAlongLine1.csv", bodyForceStats, delimiter=",")
#
## Repeat for the second line
#for i in range(numPoints):
#    # Run the placeAttachmentPoint.py script
#    print("Running placeAttachmentPoint.py with coordinates: "+str(line2[i,0])+ \
#            " "+str(line2[i,1])+" "+str(line2[i,2]))
#    os.system("python placeAttachmentPoints.py "+str(line2[i,0])+" "+ \
#            str(line2[i,1])+" "+str(line2[i,2]))
#    # Clear the case
#    os.system("pyFoamClearCase.py .")
#    # Run the case
#    os.system("./Allrun -special -serial")
#    # Calculate and store the maximum body force
#    os.system("pvpython bodyForceStatistics.py")
#    # Sanity check, if i=13, run paraFoam
#    if i == 13:
#        os.system("paraFoam")
#    elif i == 20:
#        os.system("paraFoam")
#    # Read in the data
#    data = np.genfromtxt("./results/bodyForceStats.csv", delimiter=",", \
#            skip_header=1)
#    # Store the maximum, minimum and mean body force, columns 1, 3 and 2
#    # Row 0
#    bodyForceStats[i,0] = data[1]
#    bodyForceStats[i,1] = data[3]
#    bodyForceStats[i,2] = data[2]
#
## Save the data to a csv file
#np.savetxt("bodyForceAlongLine2.csv", bodyForceStats, delimiter=",")

# Generate colormap for plotting the data with matplotlib
# import the pokepallete package
from pokemonPalette import pokePalette
# create a pokepallete object with the name of the palette
colors = pokePalette.get_pokemon_palette("charmander", 2)

# import matplotlib
import matplotlib as mpl
import matplotlib.gridspec as gridspec

# configure matplotlib
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["mathtext.fontset"] = "cm"
plt.rc('font', family='serif', size=12)

# Read in the data from the files
bodyForceStats1 = np.genfromtxt("bodyForceAlongLine1.csv", delimiter=",")
bodyForceStats2 = np.genfromtxt("bodyForceAlongLine2.csv", delimiter=",")

# Normalize the mean body force (column 2) by the maximum mean body force
maxMeanBodyForce = np.max([np.max(bodyForceStats1[:, 2]), np.max(bodyForceStats2[:, 2])])
bodyForceStats1[:, 2] = bodyForceStats1[:, 2] / maxMeanBodyForce
bodyForceStats2[:, 2] = bodyForceStats2[:, 2] / maxMeanBodyForce

# Create a custom grid layout with 2 rows and 1 column for two subplots
fig = plt.figure(figsize=(6, 6))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])

# First subplot
ax1 = plt.subplot(gs[0])
ax1.plot(line1Length, bodyForceStats1[:, 2], 'o', color='black',\
        label=r'$\~\mathbf{f}_L^{\mathrm{mean}}$', linewidth=2, linestyle='solid')
ax1.set_xticks([0, 0.077, 0.113, 0.127])
ax1.set_xticklabels(['A', 'B', 'C', 'D'])
#ax1.legend(loc='upper left', frameon=True)
ax1.tick_params(direction='in', which='both', bottom=True, top=True, left=True, right=True)
ax1.set_ylim([0.4, 1.1])

# Second subplot
ax2 = plt.subplot(gs[1])  # Share the x-axis with the first subplot
ax2.plot(line2Length, bodyForceStats2[:, 2], 'o', color='black',\
        label=r'$\mathbf{f}_L^{\mathrm{mean}}$', linewidth=2, linestyle='solid')
ax2.set_xticks([0, 0.077, 0.113, 0.127])
ax2.set_xticklabels(['H', 'G', 'F', 'E'])
#ax2.legend(loc='upper left', frameon=True)
ax2.tick_params(direction='in', which='both', bottom=True, top=True, left=True, right=True)
ax2.set_ylim([0.4, 1.1])

fig.text(0.03, 0.49, r'$||\bar{\mathbf{f}}_L||$', va='center', \
        rotation='vertical', fontsize=18)

# On each subplot create regions of color to indicate sections of the
# path. Each subplot has three regions. The first region from x=0 to x=0.077 is
# colored ed9494ff, the second region from x=0.077 to x=0.113 is colored
# 94ed94ff, and the third region from x=0.113 to x=0.127 is colored 9494edff.
# The y values are set to 0.52 and 1.08 to cover the entire y range of the
# The alpha value is set to 0.5 to make the colors transparent.
ax1.axvspan(0, 0.077, 0.04, 0.96, facecolor='#ed9494ff', alpha=0.5)
ax1.axvspan(0.077, 0.113, 0.04, 0.96, facecolor='#94ed94ff', alpha=0.5)
ax1.axvspan(0.113, 0.127, 0.04, 0.96, facecolor='#9494edff', alpha=0.5)
ax2.axvspan(0, 0.077, 0.04, 0.96, facecolor='#ed9494ff', alpha=0.5)
ax2.axvspan(0.077, 0.113, 0.04, 0.96, facecolor='#94ed94ff', alpha=0.5)
ax2.axvspan(0.113, 0.127, 0.04, 0.96, facecolor='#9494edff', alpha=0.5)

# Adjust the spacing between subplots
plt.subplots_adjust(hspace=0.3)

# Save the figure
plt.savefig("bodyForceMeanAlongLine.pdf", bbox_inches='tight', dpi=300)
