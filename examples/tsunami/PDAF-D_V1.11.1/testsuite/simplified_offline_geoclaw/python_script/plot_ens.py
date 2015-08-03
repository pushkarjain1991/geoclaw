import matplotlib.pyplot as plt
import plotmap
import matplotlib as m
import numpy as np


def plot_ens(x,y, numens):
    cdict = {
      'red'  :  ( (0.0, 0.25, .25), (0.02, .59, .59), (1., 1., 1.)),
      'green':  ( (0.0, 0.0, 0.0), (0.02, .45, .45), (1., .97, .97)),
      'blue' :  ( (0.0, 1.0, 1.0), (0.02, .75, .75), (1., 0.45, 0.45))
    }

    cm = m.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)
    #plt.figure()
    numrows = 3
    numcols = 3
    fig,axarray = plt.subplots(numrows,numcols)
    for i,a in enumerate(axarray.flatten()):
        #z = np.loadtxt("ens_0" + str(i+1) + "_ana.txt")
        z = np.loadtxt("first_ens_" + str(i+1) + ".txt")
        origin = 'lower'
        cs = a.contourf(x,y,z,cmap = cm, origin=origin, vmin=0, vmax = 100)
        #cs = a.contourf(x,y,z,cmap = plt.cm.rainbow, origin=origin, vmin=0, vmax = 100)
        a.set_title('Ensemble'+str(i+1))
        #plotmap.docontour(x,y,z)
        

    cbar = fig.add_axes([.92,0.15,0.05,0.7])
    #cbar = fig.add_axes([0.85,0.15,0.05,0.7])
    fig.colorbar(cs, cax=cbar)
    plt.show()

if __name__ == '__main__':
    x = np.linspace(-50,50,51)
    y = np.linspace(-50,50,51)
    plot_ens(x,y,9)
