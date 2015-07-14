import numpy as np

def makeinitens(z1, num_ens, outfile_ens):
    """
    Write output for ensemble members. This format is input for PDAF
    """
    
    for i in range(1,num_ens + 1):
        filename = outfile_ens + str(i) + ".txt"
        np.savetxt(filename, z1*(i-(num_ens-1.0)/2.0), fmt='%-7.5f')

