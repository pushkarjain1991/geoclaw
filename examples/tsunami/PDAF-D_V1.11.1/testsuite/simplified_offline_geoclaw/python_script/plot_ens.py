import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import plotmap
import matplotlib as m
import matplotlib.cm as cm
import numpy as np


def plot_ens(x,y, numens, initial_state = "False", analysis_state = "True"):
    numrows = 3
    numcols = 3
    fig,axarray = plt.subplots(numrows,numcols)
    origin = 'lower'
    #fig = plt.figure()
    #grid = ImageGrid(fig,111,nrows_ncols=(numrows,numcols),axes_pad=0.1)
    for i,a in enumerate(axarray.flatten()):
        if analysis_state:
            z = np.loadtxt("ens_0" + str(i+1) + "_ana.txt")
        else: 
            z = np.loadtxt("first_ens_" + str(i+1) + ".txt")
        vmax = z.max()
        vmin = z.min()

        #For radial bowl test case, use the following vmin and vmax. This will be automated in future.
        vmax = 0.6
        vmin = -0.002

        levels = np.arange(vmin, vmax, (vmax-vmin)/10)
        norm = cm.colors.Normalize(vmax=vmax,vmin=vmin)
        #cs = plotmap.docontour(x,y,z)
        cs = a.contourf(x,y,z,levels, norm = norm, cmap = plt.get_cmap('bwr'), origin=origin, vmin=vmin, vmax=vmax)
        plt.colorbar(cs,ax=a,shrink=0.9)
        a.set_title('Ensemble'+str(i+1))
        

    #cbar = fig.add_axes([.92,0.15,0.05,0.7])
    #cbar = fig.add_axes([0.85,0.15,0.05,0.7])
    #fig.colorbar(a, cax=cbar)
    if analysis_state:
        plt.suptitle("Analysis state")
    else:
        plt.suptitle("Initial state")
    
    plt.show()

if __name__ == '__main__':
    x = np.linspace(-100.0,100.0,51)
    y = np.linspace(-100.0,100.0,51)
    plot_ens(x,y,9)
