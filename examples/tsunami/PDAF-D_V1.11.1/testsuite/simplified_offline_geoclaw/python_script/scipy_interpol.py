__author__ = "Pushkar Kumar jain"

import numpy as np
from scipy import interpolate

def interpol(xx,yy,zz):
    dx = xx[0,1]-xx[0,0]
    dy = yy[1,0] - yy[0,0]
    f = interpolate.interp2d(xx,yy,zz, kind='linear')
    #xnew = np.arange(np.amin(xx) - dx/2.0, np.amax(xx) + 3*dx/2.0, dx)
    #ynew = np.arange(np.amin(yy) - dy/2.0, np.amax(yy) + 3*dy/2.0, dy)
    xnew = np.linspace(-2.5, 2.5 , 5)
    ynew = np.linspace(2.5, -2.5, 5)
    znew = f(xnew,ynew)
    return znew

if __name__=="__main__":
    oldxv = np.linspace(-2.0,2.0, 4)
    oldyv = np.linspace(2.0,-2.0, 4)
    oldxx,oldyy = np.meshgrid(oldxv,oldyv)
    oldzz = oldxx*oldxx + oldyy*oldyy
    print "input array is"
    print oldzz
    output_array = interpol(oldxx,oldyy,oldzz)
    print "output array is"
    print output_array
