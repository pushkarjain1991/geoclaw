import numpy as np
import matplotlib.pyplot as plt
import ReadAmrForLevel as ramr


def create_mask(xv,yv,zmin):
    #z = np.zeros_like(xv, dtype="int32")

    #z[xv**2 + yv**2 <= zmin*100.0+16.0] = 1

    #plt.contourf(xv,yv,z)
    #plt.show()

    #np.savetxt("masked_topo.txt", z)
    
    state = np.loadtxt("ens_" + str(j) + ".txt")
    masked_val = np.loadtxt("masked_topo.txt")
    #initial_water = np.ma.array(initial_state, mask = original_case.land==0.0)
    #initial_land = original_case.land
    initial_water = np.ma.array(initial_state, mask = masked_val==0)
    initial_land = np.ma.array(initial_state, mask = masked_val==1)
        
    plotmap.docontour(mxv,myv,initial_water, initial_land, "Ensemble # " + str(j), -2.6, 50., colorbar = True)

def create_mask2(ictype,mask_type=1):
    if mask_type == 1:
        original_case = ramr.ReadAmrForLevel("../_output_original_" + ictype +"/fort.q0000",1)
        mask_state = np.zeros_like(original_case.water,dtype="int32")
        mask_state[original_case.land==0.0] = 1

    if mask_type == 2:
        mask_state = np.zeros_like(xv, dtype="int32")
        mask_state[xv**2 + yv**2 <= zmin*100.0] = 1
        plt.contourf(xv,yv,mask_state)
        plt.show()
        
    np.savetxt("masked_topo.txt",mask_state)    
    #plotmap.docontour(mxv,myv,initial_water, initial_land, "Ensemble # " + str(j), -2.6, 50., colorbar = True)
    return mask_state


if __name__=="__main__":
    xlower = -98.0
    xupper = 98.0
    ylower = -98.0
    yupper = 98.0
    nx = 100
    ny = 100
    zmin = 80.
    ictype = "hump"

    X = np.linspace(xlower, xupper, nx)
    Y = np.linspace(ylower, yupper, ny)

    xv, yv = np.meshgrid(X,Y)

    create_mask2(ictype,mask_type=2) 

    




