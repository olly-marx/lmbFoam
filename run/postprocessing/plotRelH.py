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

# get the saveName from the command line
saveName = sys.argv[1]

# configure matplotlib
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["mathtext.fontset"] = "cm"
plt.rc('font', family='serif', size=12)

fig = plt.figure()

# get data from file ./u_rms_vs_t+saveName+.csv ignoring the first row
data = \
np.genfromtxt('/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/postprocessing/'+saveName+'.csv', delimiter=',', skip_header=1)

# plot the data time is column 0, height_i is column 3. Take the square root of the data
plt.plot(data[:,0], data[:,4], label=r'$h_E/h_0$')

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
plt.savefig('/home/ojm40/foam/ojm40-5.0/work/liquidMetalBattery/run/postprocessing/h_E-'+saveName+'.pdf',\
        bbox_inches='tight', dpi=300)

