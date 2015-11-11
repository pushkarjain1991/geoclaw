import os
import pdb
import subprocess
import shutil
import numpy as np
import make_init_ens
import pdaf_to_geoclaw
import obs
import cmdir
import maketopo
import plot_ens
import plotmap
import ReadAmrForLevel as ramr
import run_geoclaw
import verification_plot
from clawpack.geoclaw import topotools
#import myplots
# import ensemble_class
# import matplotlib.pyplot as plt
# import read_amr


def main():
    """
    Performs forecast and assimilation steps of geoclaw
    """
    amr_max_level = 1
    output_times = 12
    DA = False
    num_ens = 9
    #  dtobs = [0.0, 8.0]
    #  dtobs = [0.0, 4.0, 8.0]
    #  dtobs = [0.0, 2.0, 4.0, 6.0, 8.0]
    #  dtobs = [0.0, 2.0, 4.0, 6.0, 8.0]
    #  dtobs = [0.0, 4.0, 8.0, 12.0, 16.0]
    dtobs = np.linspace(0.0, 16.0, 5)

    #  Model Parameters
    nxpoints = 51
    mx = nxpoints-1
    nypoints = 51
    my = nypoints-1
    xlower = -100.e0
    xupper = 100.e0
    yupper = 100.e0
    ylower = -100.e0
    geoclaw_exec = "../xgeoclaw"
    topo_path = "../bowl.topotype2"

    #  DA parameters
    stddev_obs = 0.1
    dxobs = 10
    dyobs = 10
    PDAF_executable = "./PDAF_offline"

    x = np.linspace(xlower, xupper, nxpoints)
    y = np.linspace(yupper, ylower, nypoints)
    xv, yv = np.meshgrid(x, y)

    dx = (xupper-xlower)/(nxpoints-1)
    dy = (yupper-ylower)/(nypoints-1)

    x_cell = np.linspace(xlower + dx/2.0, xupper - dx/2.0, mx)
    y_cell = np.linspace(yupper - dy/2.0, ylower + dy/2.0, my)
    mxv, myv = np.meshgrid(x_cell, y_cell)

    firsttime = True

    # Create main initial data file in format of the Geoclaw input
    # Geoclaw input format
    #  ___ ___ ___
    #  ___ ___ ___
    #  ___ ___ ___
    #  ___ ___ ___
    #  ___ ___ ___
    #  ___ ___ ___
    print "Creating the initial Gaussian hump ..."
    mean_init_z = maketopo.qinit(xv, yv)
    np.savetxt("first_ensemble_main", mean_init_z)
    
    #  Original Geoclaw
    print "Running original geoclaw ..."
    cmdir.take("original_geoclaw")
    os.chdir("original_geoclaw")
    pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, "../first_ensemble_main", "hump.xyz")
    #topotools.topo1writer("hump.xyz", maketopo.qinit, xlower, xupper, ylower,
    #                      yupper, nxpoints, nypoints)
    original_timeinterval = [dtobs[0], dtobs[-1]]
    original_output_times = output_times*(np.size(dtobs)-1)
    run_geoclaw.run_geoclaw(0, dtobs=original_timeinterval, mx=mx, my=my,
                            geoclaw_input="hump.xyz",
                            topography=topo_path,
                            output_times=original_output_times,
                            max_amr=amr_max_level,
                            geoclaw_exec=geoclaw_exec)
    os.chdir("../")

    # Create ensemble members based on the mean value vector
    # Ensemble data format
    #  __ __ __ __ __ __ __ ... __
    #  __ __ __ __ __ __ __ ... __
    #  __ __ __ __ __ __ __ ... __
    #  __ __ __ __ __ __ __ ... __
    #  __ __ __ __ __ __ __ ... __
    #  __ __ __ __ __ __ __ ... __
    #  __ __ __ __ __ __ __ ... __
    #  __ __ __ __ __ __ __ ... __
    print "Creating the first ensemble members ..."
    make_init_ens.makeinitens(mean_init_z, num_ens, "first_ens_")

    for k, j in enumerate(dtobs[:-1]):

        # Write ensemble_tracker
        # For every forward run, Fortran reads the ensemble number
        # from ens_tracker.
        with open("ens_tracker", "w") as ens_tracker:
            ens_tracker.write("1")

        for i in range(1, num_ens+1):
            # Define pdaf input and output file names
            geoclaw_input = "hump_ens_" + str(i) + ".txt"
            pdaf_output = "../ens_0"+str(i)+"_ana.txt"
            subdir_name = "ens_"+str(i)
            first_ensemble = "../first_ens_" + str(i) + ".txt"
            if DA:
                pdaf_input = "../ens_" + str(i) + ".txt"
                read_geoclaw_output = "fort.q00" + str(output_times)

            #print "\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
            print "\nEnsemble unit " + str(i) + " at " + str(j) + " secs"
            #print "# # # # # # # # # # # # # # # # # # # # # # # # # # # # #"

            # Check if path exists and create subfolders to run
            # individual geoclaw
            cmdir.take(subdir_name + "_" + str(k))

            # Change to subdirectory
            os.chdir(subdir_name + "_"+str(k))

            # Prepare qinit files for geoclaw
            # Convert format of ensemble to data input format of Geoclaw qinit
            if firsttime:
                print "Input read from initial condition"
                pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, first_ensemble, \
                        geoclaw_input)
                #topotools.topo1writer(geoclaw_input, maketopo.qinit_ens,
                                      #xlower, xupper, ylower, yupper,
                                      #nxpoints, nypoints)
            else:
                if DA:
                    print "Input read from previous state"
                    pdaf_to_geoclaw.pdaf_to_geoclaw(mxv, myv, pdaf_output,
                                                    geoclaw_input)
                else:
                    np.savetxt("dummy1", np.zeros((mx, my)))
                    pdaf_to_geoclaw.pdaf_to_geoclaw(mxv, myv, "dummy1",
                                                    geoclaw_input)

            # Run Geoclaw forecast step
            # ---------------------------------------#
            # # # # # # # #         FORECAST       # # # # # # # # # #
            # ---------------------------------------#
            #  Prepare individual geoclaw input
            print "Executing forecast step for ens_number " + str(i)
            run_geoclaw.run_geoclaw(k, dtobs=dtobs, mx=mx, my=my,
                                    geoclaw_input=geoclaw_input,
                                    topography=topo_path,
                                    output_times=output_times,
                                    max_amr=amr_max_level,
                                    geoclaw_exec=geoclaw_exec)
            print "Forecast completed for ens_number " + str(i)

            # ---------------------------------------------------------------#
            # # # # # # # #         GEOCLAW OUTPUT TO PDAF INPUT I/O       # \
# # # # # # # # #
            # ---------------------------------------------------------------#
            # np.savetxt(pdaf_input, reshaped_eta)
            if DA:
                stop_case = ramr.ReadAmrForLevel(read_geoclaw_output, 1.0)
                np.savetxt(pdaf_input, stop_case.eta_with_land)

            #  Very dangerous
            if not (j == dtobs[-2]):
                shutil.copy2("fort.q0012", "../fort.q0012_" + subdir_name)
                # pdb.set_trace()
            # Go back one directory
            os.chdir("../")

        firsttime = False
        os.remove("ens_tracker")

        if DA:
            print "Data assimilation is on ...\n"
            print "Recevied observation at " + str(dtobs[k+1]) 

            # Create observation data
            print "Creating observation data ...\n"
            # mean_init_z = np.loadtxt("ens_"+str((num_ens+1)/2)+".txt")
            original_fortq_time = str((np.size(dtobs)-1) * output_times *
                                  int(dtobs[k+1])/int(dtobs[-1])).zfill(2)
            original_fortq_file = "./original_geoclaw/fort.q00"+original_fortq_time
            print "\nTrue state read is " + original_fortq_file
            original_case = ramr.ReadAmrForLevel(original_fortq_file, 1.0)

            truefield = original_case.get_water()
            observation = obs.make_obs(mx, my, dxobs, dyobs, stddev_obs,
                                       truefield)
        # Run PDAF assimilation step
        # -------------------------------------------#
        # # # # # # # #         ASSIMILATION       # # # # # # # # # #
        # -------------------------------------------#
            print "Executing assimilation step ..."
            subprocess.call(PDAF_executable, stdout=open(os.devnull,'w'),
                            stderr=subprocess.STDOUT)
            print "Assimilation completed at time " + str(dtobs[k+1]) + " secs"
            print "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n"

    # -------------------------------------------#
    # # # # # # # # # # #      POST-PROCESSING    # # # # # # # # # #
    # -------------------------------------------#
    # Verification of geoclaw restart
    verification_plot.verify_das(dtobs, (num_ens+1)/2, output_times, num_ens)

    # Plotting ensemble states
    if (num_ens != 1):
        plot_ens.plot_ens(xv, yv, num_ens, initial_state=True,
                          analysis_state=False)

    if DA:
        # Plotting observation data
        plotmap.docontour(mxv, myv, observation, stop_case.land,
                          "Observation data", -1000., 500)
        # Plotting assimilated state
        assimilated_state = np.loadtxt("state_ana.txt")
        assimilated_water = np.ma.array(assimilated_state,
                                        mask=stop_case.land == 0.0)
        np.savetxt("ha1", assimilated_water)
        # print stop_case.land
        # print stop_case.water
        assimilated_land = stop_case.land
        plotmap.docontour(mxv, myv, assimilated_water,
                          assimilated_land, "Assimilation state", -0.9, 0.9)

        # Error between assimilated state and original
        # original_fortq_time=str((np.size(dtobs)-1)*output_times*\
        # int(dtobs[k+1])/int(dtobs[-1])).zfill(2)
        original_fortq_time = str((np.size(dtobs)-1)*output_times).zfill(2)
        original_fortq_file = "./original_geoclaw/fort.q00"+original_fortq_time
        original_case = ramr.ReadAmrForLevel(original_fortq_file, 1.0)
        error_eta_water = assimilated_water - original_case.eta_with_land
        np.savetxt("ha2", error_eta_water)
        error_eta_water_percent = (assimilated_water -
                                   original_case.eta_with_land) * 100.0 / \
                                   original_case.eta_with_land
        # error_eta_water_percent = (assimilated_water - \
        # original_case.eta_with_land)*100.0/original_case.eta_with_land
        plotmap.docontour(mxv, myv, error_eta_water, original_case.land,
                          "Assimilated state error", -0.01, 0.01)
        plotmap.docontour(mxv, myv, error_eta_water_percent,
                          original_case.land, "Assimilated state error %",
                          -100.0, 100.0)

if __name__ == '__main__':
    # subprocess.call("make cleanmy")
    subprocess.call("make")
    # Construct topography
    z1 = maketopo.maketopo()

    main()
