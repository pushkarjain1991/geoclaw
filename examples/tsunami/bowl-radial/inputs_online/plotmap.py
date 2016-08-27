import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from pylab import *
from numpy import linalg as LA
import numpy.ma as ma
import myplots

def docontour(x,y,z_water, z_land, plot_title, vmin, vmax, colorbar=True, savefile=None):
    #fig = plt.figure()
    #fig.set_canvas(plt.gcf().canvas)
    #from clawpack.visclaw import colormaps, geoplot
    from clawpack.visclaw import geoplot
    
    #water
    pcolor_cmap_water = geoplot.tsunami_colormap

    #Land
    pcolor_cmap_land = geoplot.land_colors

    levels=np.linspace(vmin,vmax,100)
    ax = plt.contourf(x,y,z_land,vmin = 0.0, vmax = 100.0, cmap = plt.get_cmap(pcolor_cmap_land))
    cs = plt.contourf(x,y,z_water,levels, vmin =vmin, vmax = vmax, cmap = plt.get_cmap(pcolor_cmap_water))
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(plot_title)
    #plot_text = r"L^2: " + str(LA.norm(z_water[~z_water.mask])) + "\n" + r"L^\infty: " + str(np.max(np.abs(z_water)))
    #plt.figtext(0.4, 0.045, plot_text, color='black',weight='roman',fontsize=12,bbox={'facecolor':'white', 'boxstyle':'round'}, style='italic',ha='center', va='center')
    
    if colorbar:
        cbar = plt.colorbar(cs, shrink = 0.9,format='%.0e')
        cbar.ax.set_ylabel("WSE")

    if savefile is not None:
        print "Saving ",savefile
        plt.savefig(savefile)

    plt.close()
    #plt.show()
        

def class_contour(test_case, plot_title, vmin, vmax, savefile=None):
    #from clawpack.visclaw import colormaps, geoplot
    from clawpack.visclaw import geoplot
    
    x = test_case.mxv
    y = test_case.myv
    z_water = test_case.water
    z_land = test_case.land
   
    docontour(x,y,z_water, z_land, plot_title, vmin, vmax,savefile=savefile)

 

if __name__=='__main__':
    #plotdata.clearfigures()
    #drytol = 1.e-2
    import ReadAmrForLevel as ramr
    import myplots
   
    x = np.linspace(-98,98,50)
    y = np.linspace(-98,98,50)
    mxv,myv = np.meshgrid(x,y)
    #xobs_start = 25 
    #yobs_start = 25
    #xobs_end = 30
    #yobs_end = 30
    #nxobs = 2
    #nyobs = 2

    #obs_x = np.linspace(xobs_start,xobs_end, nxobs).astype('int32')
    #obs_y = np.linspace(yobs_start,yobs_end, nyobs).astype('int32')
    #print obs_x, obs_y
    #obs_xv, obs_yv = np.meshgrid(obs_x, obs_y)


    #obs_mat = np.ones((50,50))*-999.0

    for i in range(1,6):
        #read_geoclaw_output = "_output_original/fort.q0005"
        read_geoclaw_output = "../_output_original_hump/fort.q000" + str(i)
        original_case = ramr.ReadAmrForLevel(read_geoclaw_output, 1.0)
        original_water = original_case.water
        original_land = original_case.land
        
        assimilated_state = np.loadtxt("../_output/state_step" + str(i) +"0_ana.txt")
        assimilated_water = np.ma.array(assimilated_state,mask = original_case.land==0.0)
        assimilated_land = original_case.land
        
        print np.min(assimilated_water)
        print np.max(assimilated_water)
        docontour(mxv,myv,assimilated_water, assimilated_land, "Assimilated state", -2.6, 2.3)
        docontour(mxv,myv,original_water, original_land, "Original state", -2.6, 2.3)
        plt.show()
