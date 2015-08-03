from scipy import interpolate
import numpy as np

def test_interp(mxv, myv, reshaped_eta, x, y):
    f = interpolate.RectBivariateSpline(mxv, myv, reshaped_eta)
    #f = interpolate.interp2d(mxv, myv, reshaped_eta)
    interp_eta = f(x,y)
    return interp_eta

if __name__ == "__main__":
    
    nxpoints = 5
    mx = nxpoints-1
    nypoints = 5
    my = nypoints-1
    xlower = -50.e0
    xupper = 50.e0
    yupper = 50.e0
    ylower = -50.e0
    dx = (xupper-xlower)/(nxpoints-1)
    dy = (yupper-ylower)/(nypoints-1)
    x_cell = np.linspace(xlower + dx/2.0, xupper - dx/2.0, mx)
    y_cell = np.linspace(ylower + dy/2.0, yupper - dy/2.0, my)
    x = np.linspace(xlower,xupper,nxpoints)
    y = np.linspace(ylower,yupper,nypoints)

    mxv, myv = np.meshgrid(x_cell,y_cell)
    reshaped_eta = np.identity(mx)
    print "xcell", x_cell
    print "ycell", y_cell
    print reshaped_eta

    print "xnew",x
    print "ynew",y
    haha = test_interp(x_cell, y_cell, reshaped_eta, x, y)
    #haha = test_interp(mxv, myv, reshaped_eta, x, y)
    print "interpolated values", haha

