import numpy as np
import geoclaw_input_format as gcif


def makeinit(xv, yv, outfile):
    """
    Create main initial data file and write in geoclaw initial data format

    ARGUMENTS:
        xv, yv: Input arrays of size nx*ny
        outfile: Output file to write into

    RETURNS:
        z: WSE using f(xv,yv) for Gaussian bump. Size nx*ny
    """

    # Gaussian bump
    ze = -((xv)**2 + (yv)**2)/10.
    z = np.where(ze > -10., 40.e0 * np.exp(ze), 0.)

    # Write output for initial data in GeoCLAW input format
    gcif.geoclaw_input_format(xv, yv, z, outfile)
    return z
