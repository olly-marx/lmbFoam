#plot Ux and Bx data for Hartman validation case with matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

# for each validation case, fetch the csv file and store in a list
# there will be a list for the analytical solution and a list for each numerical solution
# first fetch the numerical solution files and store Bx and Ux data in a list
# the first file is "hartmann_mhdFoam_M20" columns 1 and 4 are Bx and Ux
# the y value is in column 8
hartmannFiles = []
hartmannFiles.append("hartmann_mhdFoam_M1.csv")
hartmannFiles.append("hartmann_mhdFoam_M20.csv")
hartmannFiles.append("hartmann_lmbFoam_M20.csv")

# Generate colormap for plotting the data with matplotlib
# import the colormap from matplotlib
from matplotlib import cm
# create a colormap with 4 colors
colors = cm.get_cmap('tab10', 4)

# configure matplotlib
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["mathtext.fontset"] = "cm"
plt.rc('font', family='serif', size=12)

# create arrays to store the data
Ux = []
y = []

# loop over the files and store the data in the arrays, skip the first line
for M in [1,5,20]:

    for solver in ["mhdFoam"]:
        file = "hartmann_" + solver + "_M" + str(M) + ".csv"

        data = np.genfromtxt(file, delimiter=',', skip_header=1)
        Ux.append(data[:,4])
        y.append(data[:,8])

        # Normalise the data, Ux by max(Ux)
        Ux[-1] = Ux[-1]/max(Ux[-1])

        # Pick a marker for each M value
        if M == 1:
            mark = 'o'
            col = colors(1)
        elif M == 5:
            mark = 'D'
            col = colors(2)
        elif M == 20:
            mark = 'x'
            col = colors(3)
        
        # plot the data, first the Ux data, the Bx data in the second loop
        # the color is set by the index of the file in the list
        # the label should be the solver used and the M value
        thisLabel = file.split("_")[1] + " (M=" + file.split("_")[2].split(".")[0].replace("M","") + ")"
        plt.plot(Ux[-1], y[-1], color=col, label=thisLabel, marker=mark, markevery=10, fillstyle='none', \
                markersize=4, linestyle='none')

    # add the analytical solution
    # the analytical solution is calculated by first creating a mesh of y values
    # then calculating the analytical solution for each y value
    # the mesh is created by using the y values from -1 to 1, and 1000 points

    # create the mesh
    yMesh = np.linspace(-1, 1, 1000)
    UxAnalytical = np.zeros(len(yMesh))

    # for each M value, calculate the analytical solution
    # first calculate delta for this case
    if M == 0:
        delta = 100000000
    else:
        delta = 1/M
    y0 = 1.0

    # calculate the analytical solution for each y value
    for i in range(len(yMesh)):
        UxAnalytical[i] = (np.cosh(M) - np.cosh(yMesh[i]/delta))/(np.cosh(M) - 1)

    # plot the analytical solution
    thisLabel = "Analytical (M=" + str(M) + ")"
    plt.plot(UxAnalytical, yMesh, color='black', linestyle='-', linewidth=1)

# add the legend
plt.legend()

# add the axis labels
plt.xlabel(r'$U_x / U_0$')
plt.ylabel(r'$y$')

# set the axis limits
plt.xlim([0,1.2])
plt.ylim([-1,1])

# save the figure as a pdf, high resolution
plt.savefig("hartmannValidationU.pdf", dpi=300, bbox_inches='tight')

# clear the figure
plt.clf()

# now do the same for the Bx data
# create arrays to store the data
Bx = []
y = []

# loop over the files and store the data in the arrays, skip the first line
for M in [1,5,20]:

    for solver in ["mhdFoam"]:
        file = "hartmann_" + solver + "_M" + str(M) + ".csv"

        data = np.genfromtxt(file, delimiter=',', skip_header=1)
        Bx.append(data[:,1])
        y.append(data[:,8])

        # Normalise the data, Bx by max(Ux)
        # Bx[-1] = Bx[-1]

        # Pick a marker for each M value
        if M == 1:
            mark = 'o'
            col = colors(1)
        elif M == 5:
            mark = 'D'
            col = colors(2)
        elif M == 20:
            mark = 'x'
            col = colors(3)
        
        # plot the data, first the Ux data, the Bx data in the second loop
        # the color is set by the index of the file in the list
        # the label should be the solver used and the M value
        thisLabel = file.split("_")[1] + " (M=" + file.split("_")[2].split(".")[0].replace("M","") + ")"
        plt.plot(Bx[-1], y[-1], color=col, label=thisLabel, marker=mark, markevery=10, fillstyle='none', \
                markersize=4, linestyle='none')

    # add the analytical solution
    # the analytical solution is calculated by first creating a mesh of y values
    # then calculating the analytical solution for each y value
    # the mesh is created by using the y values from -1 to 1, and 1000 points

    # create the mesh
    yMesh = np.linspace(-1, 1, 1000)
    BxAnalytical = np.zeros(len(yMesh))

    # for each M value, calculate the analytical solution
    # first calculate delta for this case
    delta = 1/M
    y0 = 1.0

    # calculate the analytical solution for each y value
    for i in range(len(yMesh)):
        BxAnalytical[i] = -(yMesh[i]*np.sinh(M) - np.sinh(yMesh[i]/delta))/(np.cosh(M) - 1)

    # plot the analytical solution
    plt.plot(BxAnalytical, yMesh, color='black', linestyle='-', linewidth=1)

# add the legend
plt.legend()

# add the axis labels
plt.xlabel(r'$B_x / B_0$')
plt.ylabel(r'$y$')

# set the axis limits
plt.xlim([-1,1])
plt.ylim([-1,1])

# save the figure as a pdf, high resolution
plt.savefig("hartmannValidationB.pdf", dpi=300, bbox_inches='tight')

