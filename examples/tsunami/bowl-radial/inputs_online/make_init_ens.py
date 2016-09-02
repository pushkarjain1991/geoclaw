import numpy as np
import make_obs
import chunk

def makeinitens(xv, yv, num_ens, ictype="hump"):
    """
    Write output for ensemble members. This format is input for PDAF
    """
    
    param = 1.0
    print "Gaussian hump initial condition selected"
    for i in range(1,num_ens/2+1):
        z = makeinit(xv, yv, i*param, ictype = ictype)
        np.savetxt("ens_"+str(i)+".txt", z, fmt='%12.10f')

    z = makeinit(xv, yv, 0.0, ictype = ictype)
    np.savetxt("ens_" + str(num_ens/2 + 1)+".txt", z, fmt='%12.10f')

    for i in range((num_ens/2)+2, num_ens+1):
        z = makeinit(xv, yv, i*param, ictype = ictype)
        np.savetxt("ens_"+str(i)+".txt", z, fmt='%12.10f')

def makeinitens2(xv, yv, num_ens, lower=1.5, upper=3.0, ictype = "hump"):
    """
    Write output for ensemble members. This format is input for PDAF
    """
    if ictype=='hump': 
        print "Gaussian hump initial condition selected"
    else:
        print "Planewave initial condition selected"

        #pert = lambda x: (upper-lower)*np.exp(0.8)*(x-1)/(num_ens -1) + (-5.0+lower*np.exp(0.8))
    if num_ens==1:
        pert=np.array([0])
    else:
        if ictype=='hump':
            pert = 0.01*np.random.randn(num_ens) 
        else:
            pert = np.random.randn(num_ens) 
    print "pert = " , pert
    for i in range(num_ens):
        z = makeinit(xv, yv, pert[i], ictype = ictype)
        np.savetxt("ens_"+str(i+1)+".txt", z, fmt='%12.10f')


def makeinit(xv, yv, perturb, ictype = "hump"):
    """
    Create main initial data file and write in geoclaw initial data format

    ARGUMENTS:
        xv, yv: Input arrays of size nx*ny
        outfile: Output file to write into

    RETURNS:
        z: WSE using f(xv,yv) for Gaussian bump. Size nx*ny 
    """
    # Gaussian hump
    if ictype == "hump":
        ze = -((xv)**2 + (yv)**2)/10.
        z = np.where(ze>-10., (10.e0 + perturb)*np.exp(ze), 0.)
        print "Max amplitude = ", np.max(z)

    # Planewave 
    elif ictype == "planewave":
        k = 5*np.pi/100.0
        z = 0.8*np.sin(k * (xv+10.0 + perturb))

    else:
        print "Unrecognized initial condition. Choose ictype = hump or planewave"

    return z



if __name__ == "__main__":
    nx = 100
    ny = 100
    xlower = -99.0
    xupper = 99.0
    ylower = -99.0
    yupper = 99.0
    num_ens = 1
    ictype = "hump"
    ens_range = range(num_ens)
    obs_time_list = np.linspace(20,200,10,dtype='int32')
    #ictype = "planewave"

    X = np.linspace(xlower, xupper, nx)
    Y = np.linspace(ylower, yupper, ny)
    xv,yv = np.meshgrid(X,Y)
   
    #Make initial ensemble
    #makeinitens(xv, yv, num_ens, ictype)
    makeinitens2(xv, yv, num_ens,ictype = ictype)

    if ((nx > 50) & (ny > 50)):
        chunk.chunk_write(ens_range,type1="init")
