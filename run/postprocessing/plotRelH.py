# Now plot the data
# Generate colormap for plotting the data with matplotlib
# import the pokepallete package
# from pokemonPalette import pokePalette
# # create a pokepallete object with the name of the palette
# colors = pokePalette.get_pokemon_palette("squirtle", 1)

# import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys

# get the dir to look for files from the command line
direc = sys.argv[1]

# get the saveName from the command line
saveName = sys.argv[2]

# get the variable name from the command line
varName = sys.argv[3]

# Generate colormap for plotting the data with matplotlib
# import the pokepallete package
from pokemonPalette import pokePalette

# configure matplotlib
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["mathtext.fontset"] = "cm"
plt.rc('font', family='serif', size=12)

fig = plt.figure()

# get a list of the files in the directory
import os
files = os.listdir('./'+direc)

# create a pokepallete object with the name of the palette
colors = pokePalette.get_pokemon_palette("bulbasaur", len(files))

h0 = 0.031

# for each file in the directory, plot the data
for i in range(0, len(files)):
    # get data from file ./direc/filename ignoring the first row
    data = \
    np.genfromtxt('./'+direc+'/'+files[i], delimiter=',', skip_header=2)

    # each files is named direc-value-boundarycondition.csv
    # get the value from the filename
    value = files[i].split('-')[1]

    # plot the data time is column 0, height_i is column 3. Take the square root of the data
    plt.plot(data[:,0], data[:,4]/h0, label=r'$'+varName+' = $'+value, color=colors[i])

# Add a legend
# plt.legend(fontsize=12, loc='upper right')
# Add axis labels
plt.xlabel(r'$t$ (s)')
plt.ylabel(r'$h_E/h_0$')

# Move the ticks to the inside of the axes
plt.tick_params(direction='in', which='both', bottom=True, top=True, left=True, right=True)

## Set the axis limits
#plt.xlim([-0.075, 0.075])
plt.ylim([0.02999, 0.03001])

# Save the figure
plt.savefig('./h_E-'+saveName+'.pdf',\
        bbox_inches='tight', dpi=300)

