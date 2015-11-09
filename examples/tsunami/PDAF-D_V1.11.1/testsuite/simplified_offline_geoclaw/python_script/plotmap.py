import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from pylab import *

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
    gca().set_position((0.1, .15, .6, .7))
    #plt.figtext(0.35, 0.001, "Max: " + str(np.max(z_water)),color='black',weight='roman',fontsize=12,bbox={'facecolor':'white'}, style='italic')
    plt.figtext(0.4, 0.025, "Max: " + str(np.max(np.abs(z_water))),color='black',weight='roman',fontsize=12,bbox={'facecolor':'white'}, style='italic',horizontalalignment='center')
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
    import ReadAmrForLevel as ramr
    read_geoclaw_output = "./ens_5_1/fort.q0012"
    stop_case = ramr.ReadAmrForLevel(read_geoclaw_output, 1.0)

    x = np.linspace(-98,98,50)
    y = np.linspace(-98,98,50)
    mxv,myv = np.meshgrid(x,y)
    assimilated_state = np.loadtxt("state_ana.txt")
    print assimilated_state
    assimilated_water = np.ma.array(assimilated_state,mask = stop_case.land==0.0)
    print assimilated_water
    #print stop_case.land
    #print stop_case.water
    assimilated_land = stop_case.land
    docontour(mxv,myv,assimilated_water, assimilated_land, "Final state at last assimilation", -0.9, 0.9)
