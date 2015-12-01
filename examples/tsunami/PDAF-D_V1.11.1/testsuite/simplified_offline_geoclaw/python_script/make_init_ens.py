import numpy as np

def makeinitens(z1, num_ens, outfile_ens):
    """
    Write output for ensemble members. This format is input for PDAF
    """
    
    for i in range(1,num_ens + 1):
        filename = outfile_ens + str(i) + ".txt"
        perturb_value = 0.01*(i-(num_ens + 1)/2)
        z2 = np.where(z1==0.0, 0, z1 + perturb_value)
        #z2 = z1 + *(i-(num_ens + 1.0)/2.0)
        #np.savetxt(filename, z2, fmt='%-7.5f')
        np.savetxt(filename, z2)

