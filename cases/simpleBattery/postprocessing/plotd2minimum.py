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
if len(sys.argv) > 1:
    folder = sys.argv[1]
else:
    folder = ""

# create plt figure and axes and set x axis to log scale
fig = plt.figure()
ax = fig.add_subplot(111)

if folder == "Bext":
    ax.set_xscale('log')

for bcs in ["Bottom","Sidewalls"]:
    # create a list of the files in the d2Data folder
    d2DataFiles = os.listdir("postprocessing/d2Data"+folder+"/"+bcs)
    # sort the list of files
    d2DataFiles.sort()

    print(d2DataFiles)

    # create an array to store the minimum d2 values and the parameter value for each run
    d2min = np.zeros((len(d2DataFiles),2))

    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams["mathtext.fontset"] = "cm"
    plt.rc('font', family='serif', size=12)

    # loop over the files in the d2Data folder
    for file in d2DataFiles:
        # import data from the file
        data = np.genfromtxt("postprocessing/d2Data"+folder+"/"+bcs+"/"+file, delimiter=',', skip_header=1)
        # find the minimum value of d2 store it in the array
        d2min[d2DataFiles.index(file),0] = np.min(data[:,1])
        # now find the parameter value and store it in the array 
        # the parameter value is in the filename and is the first number
        d2min[d2DataFiles.index(file),1] = float(file.split("_")[1][:-4])

    #sort the array by the parameter value
    d2min = d2min[d2min[:,1].argsort()]
    # plot the data on the same set of axes
    ax.plot(d2min[:,1],d2min[:,0],label="$d_2^{min}$ "+bcs, marker='o', \
            linewidth=1, markersize=1.5)


# add a legend
plt.legend()
# add axis labels
if folder == "Ratio":
    plt.xlabel("Ratio X length to Y length")
elif folder == "Bext":
    plt.xlabel("$B_{ext}$ (T)")
plt.ylabel("$d_2$ (m)")
# save the figure
plt.savefig("results/d2GlobalMinimum_"+folder+"Study.pdf",bbox_inches='tight')
