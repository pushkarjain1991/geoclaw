"""
Create topo and dtopo files needed for this example:
    etopo10min120W60W60S0S.asc        download from GeoClaw topo repository
    dtopo_usgs100227.tt3              create using Okada model 
Prior to Clawpack 5.2.1, the fault parameters we specified in a .cfg file,
but now they are explicit below.
    
Call functions with makeplots==True to create plots of topo, slip, and dtopo.
"""

import os
import sys

import clawpack.clawutil.data
import matplotlib.pyplot as plt

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')

def get_topo(makeplots=False):
    """
    Retrieve the topo file from the GeoClaw repository.
    """
    from clawpack.geoclaw import topotools
    topo_fname = 'etopo10min120W60W60S0S.asc'
    url = 'http://www.geoclaw.org/topo/etopo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir, 
            file_name=topo_fname, verbose=True)

    if makeplots:
        from matplotlib import pyplot as plt
        topo = topotools.Topography(os.path.join(scratch_dir,topo_fname), topo_type=2)
        topo.plot()
        fname = os.path.splitext(topo_fname)[0] + '.png'
        plt.savefig(fname)
        print "Created ",fname


    
def make_dtopo(pert, filename, makeplots=False):
    """
    Create dtopo data file for deformation of sea floor due to earthquake.
    Uses the Okada model with fault parameters and mesh specified below.
    """
    from clawpack.geoclaw import dtopotools
    import numpy

    dtopo_fname = os.path.join(scratch_dir, "dtopo_usgs100227_bias.tt3_"+filename)
    #dtopo_fname = os.path.join(scratch_dir, "dtopo_usgs100227.tt3_"+filename)

    # Specify subfault parameters for this simple fault model consisting
    # of a single subfault:

    usgs_subfault = dtopotools.SubFault()
    #usgs_subfault.strike = 16. + pert[0]
    #usgs_subfault.length = 450.e3 + pert[1]
    #usgs_subfault.width = 100.e3 + pert[2]
    #usgs_subfault.depth = 35.e3 + pert[3]
    #usgs_subfault.slip = 15. + pert[4]
    #usgs_subfault.rake = 104. + pert[5]
    #usgs_subfault.dip = 14. + pert[6]
    #usgs_subfault.longitude = -72.668 + pert[7]
    #usgs_subfault.latitude = -35.826 + pert[8]
    #usgs_subfault.coordinate_specification = "top center"
    
    usgs_subfault.strike = 16.
    usgs_subfault.length = 450.e3
    usgs_subfault.width = 100.e3
    usgs_subfault.depth = 35.e3
    usgs_subfault.slip = 15. + pert
    usgs_subfault.rake = 104.
    usgs_subfault.dip = 14.
    usgs_subfault.longitude = -72.668
    usgs_subfault.latitude = -35.826
    usgs_subfault.coordinate_specification = "top center"

    fault = dtopotools.Fault()
    fault.subfaults = [usgs_subfault]

    print "Mw = ",fault.Mw()
    if(fault.Mw() == 0.0):
        sys.exist("No tsunami formed")
        

    if os.path.exists(dtopo_fname):
        print "*** Not regenerating dtopo file (already exists): %s" \
                    % dtopo_fname
    else:
        print "Using Okada model to create dtopo file"

        x = numpy.linspace(-77, -67, 100)
        y = numpy.linspace(-40, -30, 100)
        times = [1.]

        fault.create_dtopography(x,y,times)
        dtopo = fault.dtopo
        dtopo.write(dtopo_fname, dtopo_type=3)


    if makeplots:
        from matplotlib import pyplot as plt
        if fault.dtopo is None:
            # read in the pre-existing file:
            print "Reading in dtopo file..."
            dtopo = dtopotools.DTopography()
            dtopo.read(dtopo_fname, dtopo_type=3)
            x = dtopo.x
            y = dtopo.y
        plt.figure(figsize=(12,7))
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122)
        fault.plot_subfaults(axes=ax1,slip_color=True)
        ax1.set_xlim(x.min(),x.max())
        ax1.set_ylim(y.min(),y.max())
        dtopo.plot_dZ_colors(1.,axes=ax2)
        fname = os.path.splitext(os.path.split(dtopo_fname)[-1])[0] + '.pdf'
        plt.savefig(fname)
        print "Created ",fname


if __name__=='__main__':
    import numpy as np
    import sys

    if (len(sys.argv) == 1):
        sys.stdout.write('No command line dim_ens given. Using 120')
        dim_ens = 120
    else:
        dim_ens = int(sys.argv[1])
        sys.stdout.write('Using dim_ens ' + str(dim_ens))


    get_topo(False)
    #TWIN
    #var_strike = 1.6
    #var_length = 45.e3
    #var_width = 10.e3
    #var_depth = 35.e2
    #var_slip = 1.5
    #var_rake = 10.4
    #var_dip = 1.4
    #var_longitude = 2.
    #var_latitude = 2.
    
    #var_strike = 4.0
    #var_length = 150.e3
    #var_width = 50.e3
    #var_depth = 10.e3
    #var_slip = 5.0
    #var_rake = 20.0
    #var_dip = 4.0
    #var_longitude = 3.
    #var_latitude = 3.
    
    #usgs_subfault.strike = 16. + pert[0]
    #usgs_subfault.length = 450.e3 + pert[1]
    #usgs_subfault.width = 100.e3 + pert[2]
    #usgs_subfault.depth = 35.e3 + pert[3]
    #usgs_subfault.slip = 15. + pert[4]
    #usgs_subfault.rake = 104. + pert[5]
    #usgs_subfault.dip = 14. + pert[6]
    #usgs_subfault.longitude = -72.668 + pert[7]
    #usgs_subfault.latitude = -35.826 + pert[8]
    
    mu = 20.0
    ##mu = np.array([2., 50.e3, 10.e3, 4.e3, 2., 5., 2., -3., 3.])
    #sigma = np.diagflat([var_strike, var_length, var_width, var_depth, var_slip, var_rake, var_dip, var_longitude, var_latitude])
    sigma = 4.0

    np.random.seed(45)
    #s = np.random.normal(mu,sigma, dim_ens)

    #with open('pert_values.txt','w') as f1:
    #    for ens_num in range(dim_ens):
    #        filename = "ens_" + str(ens_num)
    #        pert = np.random.multivariate_normal(mu,sigma)
    #        make_dtopo(pert, filename, makeplots=False)
    #        print pert

    perts = np.random.normal(mu, sigma, dim_ens)
    for ens_num, pert in enumerate(perts):
        filename = "ens_" + str(ens_num)
        make_dtopo(pert, filename, makeplots=False)
        print pert


    # Diagnosis
    print "Sample_mean - mu: ", abs(mu - np.mean(perts))

    print "sigma - sample_std: ",abs(sigma - np.std(perts, ddof=1))

    #Plot histogram
    count, bins, ignored = plt.hist(perts, 30, normed=True)
    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ),linewidth=2, color='r')
    plt.xlabel("Sample mean")
    plt.title('$\mu=10$: Sample size = ' + str(dim_ens))
    plt.show()

