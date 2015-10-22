import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def docontour(x,y,z_water, z_land, plot_title, vmin, vmax):
    #from clawpack.visclaw import colormaps, geoplot
    from clawpack.visclaw import geoplot
    
    #water
    pcolor_cmap_water = geoplot.tsunami_colormap

    #Land
    pcolor_cmap_land = geoplot.land_colors

    levels=np.linspace(vmin,vmax,100)
    ax = plt.contourf(x,y,z_land,cmap = plt.get_cmap(pcolor_cmap_land))
    cs = plt.contourf(x,y,z_water,levels, vmin =vmin, vmax = vmax, cmap = plt.get_cmap(pcolor_cmap_water))
    plt.xlabel('X')
    plt.ylabel('Y')
    cbar = plt.colorbar(cs, shrink = 0.9)
    cbar.ax.set_ylabel("WSE")
    plt.title(plot_title)
    #cbar.add_lines
    plt.show()
    return cs



def class_contour(test_case, plot_title, vmin, vmax):
    #from clawpack.visclaw import colormaps, geoplot
    from clawpack.visclaw import geoplot
    
    x = test_case.mxv
    y = test_case.myv
    z_water = test_case.water
    z_land = test_case.land
   
    docontour(x,y,z_water, z_land, plot_title, vmin, vmax)

 

if __name__=='__main__':
    #plotdata.clearfigures()
    #drytol = 1.e-2

    #Water
    #z = np.loadtxt("state_ana.txt")
    z = np.loadtxt("eta_with_land_0.txt", unpack=True, usecols=[0])
    #z = np.loadtxt("eta_after_interp_6.0.txt", unpack=True, usecols=[0])
    print z
    np.reshape(z, (50,50))
    print z
    #z.reshape((-1,50))
    print np.shape(z)
    print z[1:5]
    if np.shape(z) == (50*50,):
        np.reshape(z, (50,50))
        print z
        x = np.linspace(-98,98,50)
        y = np.linspace(-98,98,50)
        xv,yv = np.meshgrid(x,y)
    else:
        x = np.linspace(-100,100,51)
        y = np.linspace(-100,100,51)
        xv,yv = np.meshgrid(x,y)        
        np.reshape(z, (51,51))
    fig = plt.figure()
    docontour(x,y,z)
    plt.savefig("eta_before_interp_0")
    plt.show()

