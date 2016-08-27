
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw.topotools import Topography
from numpy import *

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 201
    nypoints = 201
    xlower = -100.e0
    xupper = 100.e0
    yupper = 100.e0
    ylower = -100.e0
    outfile= "bowl.topotype2"     

    topography = Topography(topo_func=topo)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=2, Z_format="%22.15e")

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 51
    nypoints = 51
    xlower = -100.e0
    xupper = 100.e0
    yupper = 100.e0
    ylower = -100.e0
    outfile= "hump.xyz"     

    topography = Topography(topo_func=qinit)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=1)

def topo(x,y):
    """
    Parabolic bowl
    """
    # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
    zmin = 80.
    z = 1.e-2*(x**2 + y**2) - zmin
    return z


def qinit(x,y):
    """
    Gaussian hump:
    """
    from numpy import where
    ze = -((x+0e0)**2 + (y+0e0)**2)/10.
    #z = where(ze>-10., 40.e0*exp(ze), 0.)
    z = where(ze>-10., 5.e0*exp(ze), 0.)
    return z

def makeinit_planewave(x,y):
    k = 5*pi/100.0
    z = 0.8*sin(k * (x+10.0))
    #return np.transpose(z)
    return z

if __name__=='__main__':
    maketopo()
    makeqinit()
    topotools.topo1writer("planewave.xyz",makeinit_planewave, -100.0, 100.0, -100.0, 100.0, 50,50)
