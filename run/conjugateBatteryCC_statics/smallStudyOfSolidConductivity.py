# Run ./Allrun for a number of different solid conductivities
# for each case, then run solidSliceCurrent.py script to process the results

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Set up the conductivity values to as a ratio of the Pb-Bi conductivity
# which is 781000 S/m

solidConductivity = [85.0]

# Loop over the different solid conductivities
for i in range(len(solidConductivity)):
    # sed command to change the conductivity in the 0/solid/solidConductivity file
    os.system("sed -i 's/internalField   .*/internalField   " + \
            str(solidConductivity[i]*781000) + ";/g' 0/solid/solidConductivity")
    os.system("sed -i 's/internalField   .*/internalField   " + \
            str(solidConductivity[i]*781000) + ";/g' save/solid/solidConductivity")
    # Run the case
    print("Running case with solid conductivity = " + str(solidConductivity[i]))
    os.system("rm -rf 0.*")
    os.system("conjugateLmbFoam > log")
    # Run the solidSliceCurrent.py script to process the results
    print("Processing results for case with solid conductivity = " + \
            str(solidConductivity[i]))
    os.system("pvpython solidSliceCurrent.py "+str(solidConductivity[i]))

# Set up the conductivity values to plot
solidConductivity = [0.1, 1.0, 10.0, 100.0]

# Generate colormap for plotting the data with matplotlib
# import the pokepallete package
from pokemonPalette import pokePalette
# create a pokepallete object with the name of the palette
colors = pokePalette.get_pokemon_palette("squirtle", len(solidConductivity))

# import matplotlib
import matplotlib as mpl

# configure matplotlib
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["mathtext.fontset"] = "cm"
plt.rc('font', family='serif', size=12)

fig = plt.figure()

# Loop over the different solid conductivities
for i in range(len(solidConductivity)):
    # Read in the data
    data = \
    np.genfromtxt("./results/solidCurrentAlongBase_"+str(solidConductivity[i])+".csv", \
            delimiter=",", skip_header=1)
    # Plot the data, the x axis is column 5, the y axis is the magnitude of the
    # columns 1, 2, 3 normalised by the maximum value
    # Label should be formatted to one significant figure
    l = r'$\~{\sigma}_{cc} = %.1f$' % solidConductivity[i]
    norm = np.max(np.sqrt(data[:,1]**2 + data[:,2]**2 + data[:,3]**2))
    plt.plot(data[:,4], np.sqrt(data[:,1]**2 + data[:,2]**2 + data[:,3]**2)/norm, \
            label=l, color=colors[i], linewidth=2)

# Add a legend
plt.legend(fontsize=12, loc='upper right')
# Add axis labels
plt.xlabel(r'$y$ (m)')
plt.ylabel(r'$|\mathbf{J}|$ (A/m$^2$)')

# Move the ticks to the inside of the axes
plt.tick_params(direction='in', which='both', bottom=True, top=True, left=True, right=True)

# Set the axis limits
plt.xlim([-0.075, 0.075])
plt.ylim([-0.02, 1.02])

# Save the figure
plt.savefig("solidCurrentAlongBase.pdf", bbox_inches='tight', dpi=300)
