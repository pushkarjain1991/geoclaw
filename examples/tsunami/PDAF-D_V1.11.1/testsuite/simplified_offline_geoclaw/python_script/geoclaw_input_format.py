import numpy as np


def geoclaw_input_format(xv, yv, z1, outfile):
    """Write output for initial data. This format
    is qinit format for GeoCLAW"""

    dstacked = np.dstack((xv, yv, z1))
    with file(outfile, 'w') as outfile:
        for slice_2d in dstacked:
            # np.savetxt(outfile,slice_2d, fmt='%-7.5f')
            np.savetxt(outfile, slice_2d, fmt='%4.1f %4.1f %7.5f')
