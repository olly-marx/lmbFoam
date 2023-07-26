#plot the Energy vs Y and the Height vs Y for the LMB simulations
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

# Get the data from the file in the directory
# The file is the energy data ./EvY.csv
# and the file is the height data ./HyY.csv
# The data is in the format:
# time (s), avg. energy (J), x (m), y (m), z (m)
# time (s), avg. height (m), x (m), y (m), z (m)

# Get the data from the files
# Get the energy data
dataE = np.genfromtxt('./EvY.csv', delimiter=',', skip_header=1)
# Get the height data
dataH = np.genfromtxt('./HvY.csv', delimiter=',', skip_header=1)

# Order the data by y
dataE = dataE[dataE[:,3].argsort()]
dataH = dataH[dataH[:,3].argsort()]

# Convert the y data to mm
dataE[:,3] = dataE[:,3]*1000
dataH[:,3] = dataH[:,3]*1000

# Plot the data
fig, ax1 = plt.subplots()
# Set the second y-axis label
ax2 = ax1.twinx()
# Plot the energy data
ax1.plot(dataE[:,3], dataE[:,1], color=colors[0], label='Avg. \"Energy\"')
# Plot the height data
ax1.plot(dataH[:,3], dataH[:,1], color=colors[1], label='Avg. Height')

# Set the y-axis label
ax1.set_ylabel(r'$|u|^2$ (m$^2$/s$^2$)')
ax2.set_ylabel('Height (mm)')
# Set the x-axis label
ax1.set_xlabel('y (mm)')

# Save the figure
plt.savefig('EHvY.pdf', bbox_inches='tight', dpi=300)
plt.savefig('EHvY.png', bbox_inches='tight', dpi=300)
plt.close()
