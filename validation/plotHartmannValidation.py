#plot Ux and Bx data for Hartman validation case with matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

# for each validation case, fetch the csv file and store in a list
# there will be a list for the analytical solution and a list for each numerical solution
# first fetch the numerical solution files and store Ux data in a list
# the y value is in column 8

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

# create arrays to store the data
Ux = []
y = []

# create 2D array to store the L1 and L2 norms for each solver and number of cells
errorDataLmb = []
errorDataMhd = []

# loop over the files and store the data in the arrays, skip the first line
for Ha in [1, 5, 20]:

    # add the analytical solution
    # the analytical solution is calculated by first creating a mesh of y values
    # then calculating the analytical solution for each y value
    # the mesh is created by using the y values from -1 to 1, and 1000 points

    # create the mesh
    yMesh = np.linspace(-1, 1, 1000)
    UxAnalytical = np.zeros(len(yMesh))

    # for each Ha value, calculate the analytical solution
    # first calculate delta for this case
    if Ha == 0:
        delta = 100000000
    else:
        delta = 1/Ha
    y0 = 1.0

    # calculate the analytical solution for each y value
    for i in range(len(yMesh)):
        UxAnalytical[i] = (np.cosh(Ha) - np.cosh(yMesh[i]/delta))/(np.cosh(Ha) - 1)

    # plot the analytical solution
    thisLabel = "$Ha$=" + str(Ha)
    plt.plot(UxAnalytical, yMesh, color='black', linestyle='-', linewidth=1)

    # place a label on the plot to the right of the analytical solution depending on the Ha value
    if Ha == 1:
        plt.text(0.6, 0.5, thisLabel, horizontalalignment='left', verticalalignment='center')
    elif Ha == 5:
        plt.text(0.82, 0.72, thisLabel, horizontalalignment='left', verticalalignment='center')
    elif Ha == 20:
        plt.text(0.97, 0.9, thisLabel, horizontalalignment='left', verticalalignment='center')


    for solver in ["mhdFoam", "lmbFoam"]:
        for NCells in [20, 40, 80, 160]:
            file = "hartmann_" + solver + "_Ha" + str(Ha) + "_NCells" + str(NCells) + ".csv"

            data = np.genfromtxt(file, delimiter=',', skip_header=1)
            Ux.append(data[:,0])
            y.append(data[:,4])

            # Normalise the data, Ux by max(Ux)
            Ux[-1] = Ux[-1]/max(Ux[-1])

            # Pick a marker for each M value
            if Ha == 1:
                # circle
                mark = 'o'
                colorindex = 0
            elif Ha == 5:
                # diamond
                mark = 'D'
                colorindex = 2
            elif Ha == 20:
                # cross
                mark = 'x'
                colorindex = 4
            
            # Pick a color for each solver
            if solver == "mhdFoam":
                col = colors[colorindex]
            elif solver == "lmbFoam":
                col = colors[colorindex+1]

            # Define functions to calculate the L1 and L2 norms
            def l1_norm(error, NCells):
                return np.sum(np.abs(error))/NCells

            def l2_norm(error, NCells):
                return np.sqrt(np.sum(error**2)/NCells)

            # Calculate the L2 norm of the error between the analytical solution and the numerical solution
            # first interpolate the analytical solution to the y values of the numerical solution
            UxAnalyticalInterp = np.interp(y[-1], yMesh, UxAnalytical)
            # then calculate the error
            error = Ux[-1] - UxAnalyticalInterp

            # Calculate the L2 norm of the error between the analytical solution and the numerical solution
            L1norm = l1_norm(error, NCells)
            L2norm = l2_norm(error, NCells)

            # Print Solver, NCells, Ha, L1 norm and L2 norm
            print(solver + ", Res:" + str(NCells) + ", Ha:" + str(Ha) \
                    + "\nL1 norm:" + str(L1norm) + "\nL2 norm:" + str(L2norm))


            # append the an array with the Ha value, number of cells, L1 norm and L2 norm
            # to the 2D array
            if solver == "mhdFoam":
                errorDataMhd.append([Ha, NCells, L1norm, L2norm])
            elif solver == "lmbFoam":
                errorDataLmb.append([Ha, NCells, L1norm, L2norm])

            # plot the data, first the Ux data only for NCells=40
            if NCells == 40:
                thisLabel = file.split("_")[1] + " ($Ha=$" + \
                        file.split("_")[2].split(".")[0].replace("Ha","") + ")"
                plt.plot(Ux[-1], y[-1], color=col, label=thisLabel, marker=mark, \
                        fillstyle='none', markersize=3, linestyle='none')

# add the legend, where the analytical solution not included
plt.legend()

# add the axis labels
plt.xlabel(r'$u_x / u_0$')
plt.ylabel(r'$y$')

# set the axis limits
plt.xlim([0,1.2])
plt.ylim([-1,1])

# show the figure
#plt.show()

# save the figure as a pdf, high resolution
plt.savefig("hartmannValidationU.pdf", dpi=300, bbox_inches='tight')

# Clear figure
plt.clf()

# Now, on the same set of axes plot the L1 and L2 norm for each solver and with NCells
# as the x axis. Do so for each Ha value, save the plot as a pdf

for Ha in [20]:
    # create a list of the L1 and L2 norms for this Ha value, for each solver
    L1normsMhd = []
    L2normsMhd = []
    L1normsLmb = []
    L2normsLmb = []
    NCells = []
    for i in range(len(errorDataMhd)):
        if errorDataMhd[i][0] == Ha:
            L1normsMhd.append(errorDataMhd[i][2])
            L2normsMhd.append(errorDataMhd[i][3])
            NCells.append(errorDataMhd[i][1])
    for i in range(len(errorDataLmb)):
        if errorDataLmb[i][0] == Ha:
            L1normsLmb.append(errorDataLmb[i][2])
            L2normsLmb.append(errorDataLmb[i][3])

    # Make the x axis a numpy array
    NCells = np.array(NCells)

    # Now calculate reference lines for the L1 and L2 norms
    # The slopes should be -1 and -2 in the log-log plot
    # The reference values are set such that the reference lines are not obscured
    # by the data points
    FirstOrderSlope = -1
    SecondOrderSlope = -2

    # calculate the reference values
    FirstOrderRef = 0.001
    SecondOrderRef = 0.01
    
    # calculate the reference lines such that the value for nCells = 20 is FirstOrderRef
    # and the value for nCells = 20 is SecondOrderRef
    FirstOrderLine = FirstOrderRef*(NCells/20)**FirstOrderSlope
    SecondOrderLine = SecondOrderRef*(NCells/20)**SecondOrderSlope
    

    # plot the reference lines
    plt.plot(NCells, FirstOrderLine, color='k', linestyle=':')
    plt.plot(NCells, SecondOrderLine, color='k', linestyle=':')

    # place a label on the plot next to the reference lines
    plt.text(0.09, 0.66, r'$\mathcal{O}(2^{nd})$', fontsize=12, transform=plt.gca().transAxes)
    plt.text(0.09, 0.33, r'$\mathcal{O}(1^{st})$', fontsize=12, transform=plt.gca().transAxes)

    # scale the x and y axes logarithmically
    plt.xscale('log')
    plt.yscale('log')

    # plot the L1 norm
    plt.plot(NCells, L1normsMhd, color=colors[0], linestyle='-', marker='o', \
            fillstyle='none', markersize=3, label=r'$L_1$ (mhdFoam)')
    plt.plot(NCells, L1normsLmb, color=colors[1], linestyle='-', marker='o', \
            fillstyle='none', markersize=3, label=r'$L_1$ (lmbFoam)')

    # plot the L2 norm
    plt.plot(NCells, L2normsMhd, color=colors[2], linestyle='--', marker='o', \
            fillstyle='none', markersize=3, label=r'$L_2$ (mhdFoam)')
    plt.plot(NCells, L2normsLmb, color=colors[3], linestyle='--', marker='o', \
            fillstyle='none', markersize=3, label=r'$L_2$ (lmbFoam)')

    # add the legend
    plt.legend()

    # add the axis labels
    plt.xlabel(r'$N_{cells}$')
    plt.ylabel(r'$L_1$ and $L_2$ norms')

    # set the axis limits
    plt.xlim([10,300])
    plt.ylim([1e-4,1e-1])

    # save the figure as a pdf, high resolution
    plt.savefig("hartmannValidationL1L2_Ha" + str(Ha) + ".pdf", dpi=300, bbox_inches='tight')

    # Clear figure
    plt.clf()

