import matplotlib.pyplot as plt
import numpy as np

def docontour(x,y,z):
    from clawpack.visclaw import colormaps, geoplot
    
    #plotData = ClawPlotData()
    
    #water
    pcolor_cmap = geoplot.tsunami_colormap

    #Land
    pcolor_cmap = geoplot.land_colors

    origin = 'lower'
    #cs = plt.contourf(x,y,z,10,cmap = pcolor_cmap,origin=origin)
    cs = plt.contourf(x,y,z,10,origin=origin, vmin = -0.05, vmax = 0.1)
    #cs = plt.contourf(x,y,z,10,cmap = plt.cm.rainbow,origin=origin, vmin=-50, vmax = 90)
    #cs2 = plt.contour(cs, levels = cs.levels[::2], colors='r', origin=origin, hold='on')
    plt.xlabel('X')
    plt.ylabel('Y')
    cbar = plt.colorbar(cs, shrink = 0.9)
    cbar.ax.set_ylabel("WSE")
    #cbar.add_lines
    plt.show()
    return cs


if __name__=='__main__':
    #plotdata.clearfigures()
    #drytol = 1.e-2

    #Water
    #z = np.loadtxt("state_ana.txt")
    z = np.loadtxt("eta_with_land_0.txt", unpack=True, usecols=[0])
    print z
    np.reshape(z, (50,50))
    print z
    #z.reshape((-1,50))
    print np.shape(z)
    print z[1:5]
    if np.shape(z) == (50*50,):
        np.reshape(z, (50,50))
        print z
        x = np.linspace(-50,50,50)
        y = np.linspace(-50,50,50)
        xv,yv = np.meshgrid(x,y)
    else:
        x = np.linspace(-50,50,51)
        y = np.linspace(-50,50,51)
        xv,yv = np.meshgrid(x,y)        
        np.reshape(z, (51,51))
    fig = plt.figure()
    docontour(x,y,z)
    plt.savefig("eta_before_interp_0")
    plt.show()

