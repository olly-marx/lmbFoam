#plot d2 data with matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

# for each d2 run csv found in postprocessing/d2Data folder plot the data on the
# same set of axes with matplotlib
# start by looping over files in d2Data folder and importing data

# Take an armument from the command line to specify the folder to plot
# if no argument is given, use the default folder
import sys
if len(sys.argv) > 2:
    folder = sys.argv[1]
    bcs = sys.argv[2]
else:
    folder = ""

# create a list of the files in the d2Data folder
d2DataFiles = os.listdir("postprocessing/d2Data"+folder+"/"+bcs)
# sort the list of files
d2DataFiles.sort()

print(d2DataFiles)

# Generate colormap for plotting the data with matplotlib
# import the colormap from matplotlib
from matplotlib import cm
# create a colormap with the same number of colors as there are files
colors = cm.get_cmap('viridis', len(d2DataFiles))

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["mathtext.fontset"] = "cm"
plt.rc('font', family='serif', size=12)

# loop over the files in the d2Data folder
for file in d2DataFiles:
    # import data from the file
    data = np.genfromtxt("postprocessing/d2Data"+folder+"/"+bcs+"/"+file, delimiter=',', skip_header=1)
    # sort the data by time using numpy sort algorithm
    data = data[data[:,0].argsort()]
    # choose a color for the plot
    color = colors(d2DataFiles.index(file))
    # Get the index of the first "t" in the file name to use as the label
    pos1 = file.find("t")
    # Cut the file name after this index
    file = file[pos1+1:]
    # now cut the extension from the end of the file name
    file = file[:-4]
    # Finally replace the underscores with spaces
    file = file.replace("_", " ")
    # plot the data
    plt.plot(data[:,0], data[:,1], label=file, color=color, linewidth=1.5, marker='o', markersize=3)

# add a legend with smaller font
plt.legend(fontsize=10)
# add axis labels
plt.xlabel("time (s)")
plt.ylabel("$d_2$ (m)")
# save the figure
plt.savefig("results/d2AllRuns_"+folder+bcs+"Study.png",bbox_inches='tight')
