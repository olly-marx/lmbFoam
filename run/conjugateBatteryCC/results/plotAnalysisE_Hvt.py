#plot the Energy and Height vs time for the LMB simulations
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

# Generate colormap for plotting the data with matplotlib
# import the pokepallete package
from pokemonPalette import pokePalette
# create a pokepallete object with the name of the palette
colors = pokePalette.get_pokemon_palette("bulbasaur", 6)

# configure matplotlib
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["mathtext.fontset"] = "cm"
plt.rc('font', family='serif', size=12)

# Get the data from the two files in the directory
# The first file is the energy data ./EvT.csv
# The second file is the height data ./Hvt.csv
# The data is in the format:
# time (s), avg. energy (J), max energy (J)
# time (s), avg. height (m), max height (m)

# Get the data from the files
# Get the energy data
dataE = np.genfromtxt('./EvT.csv', delimiter=',', skip_header=1)
# Get the height data
dataH = np.genfromtxt('./HvT.csv', delimiter=',', skip_header=1)

# Plot the data on a single plot with two y-axes
fig, ax1 = plt.subplots()
# Plot the energy data
ax1.plot(dataE[:,0], dataE[:,1], color=colors[0], label='Avg. \"Energy\"')
#ax1.plot(dataE[:,0], dataE[:,2], color=colors[1], label='Max Energy')
# Set the y-axis label
ax1.set_ylabel(r'$|u|^2$ (m$^2$/s$^2$)')
# Set the x-axis label
ax1.set_xlabel('Time (s)')

# Create a second y-axis
ax2 = ax1.twinx()
# Plot the height data in mm
ax2.plot(dataH[:,0], dataH[:,1]*1000, color=colors[2], label='Avg. Height')

# Set the y-axis label
ax2.set_ylabel(r'$H_E$ (mm)')

# Add one legend for all the lines
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='upper left')

# Save the figure
plt.savefig('EvT_Hvt.pdf', bbox_inches='tight', dpi=300)
plt.savefig('EvT_Hvt.png', bbox_inches='tight', dpi=300)
plt.close()
