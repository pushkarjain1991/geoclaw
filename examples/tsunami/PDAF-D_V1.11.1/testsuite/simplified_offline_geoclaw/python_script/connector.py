import os
import pdb
import sys
import subprocess
import shutil
import numpy as np
import os
import sys
import make_init_ens
import pdaf_to_geoclaw
import obs
import cmdir
import ensemble_class 
import maketopo
import plot_ens
import plotmap
import matplotlib.pyplot as plt
import read_amr
import ReadAmrForLevel as ramr
import run_geoclaw
import verification_plot
from clawpack.geoclaw import topotools


def main():
    """
    Performs forecast and assimilation steps of geoclaw
    """
    amr_max_level = 1
    output_times = 12
    DA = False
    num_ens = 9
    dtobs = [0.0, 2.0, 4.0, 6.0, 8.0]
    #dtobs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    #dtobs = [0.0,8.0]
    #dtobs = [0.0,4.0,8.0]
    
    
    #Model Parameters
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
    
    #DA parameters
    stddev_obs = 0.1
    dxobs = 10
    dyobs = 10
    PDAF_executable = "./PDAF_offline"
   
    x = np.linspace(xlower,xupper,nxpoints)
    y = np.linspace(yupper,ylower,nypoints)
    xv,yv = np.meshgrid(x,y)
    
    dx = (xupper-xlower)/(nxpoints-1)
    dy = (yupper-ylower)/(nypoints-1)
    
    x_cell = np.linspace(xlower + dx/2.0, xupper - dx/2.0, mx)
    y_cell = np.linspace(yupper - dy/2.0, ylower + dy/2.0,my)
    mxv, myv = np.meshgrid(x_cell,y_cell)


    firsttime = True
    
    
    #Original Geoclaw
    print "Running original geoclaw ..."
    cmdir.take("original_geoclaw")
    os.chdir("original_geoclaw")
    #mean_init_z = maketopo.qinit(xv/2.0, yv/2.0)
    #make_init_ens.makeinitens(mean_init_z,1, "first_ens_")
    #pdaf_to_geoclaw.pdaf_to_geoclaw(xv/2.0, yv/2.0, infile="first_ens_1.txt", outfile="hump_ens_1.txt")
    topotools.topo1writer("hump.xyz",maketopo.qinit, xlower, xupper, ylower, yupper, nxpoints,nypoints)
    original_timeinterval = [dtobs[0],dtobs[-1]]
    original_output_times = output_times*(np.size(dtobs)-1)
    run_geoclaw.run_geoclaw(0, dtobs= original_timeinterval, mx = mx, my = my, geoclaw_input =  "hump.xyz", topography = topo_path, output_times = original_output_times, max_amr = amr_max_level, geoclaw_exec = geoclaw_exec)
    os.chdir("../")
    
    
    #Create main initial data file in format of the Geoclaw input
    #Create hump.xyz
    #Geoclaw input format
    # ___ ___ ___
    # ___ ___ ___
    # ___ ___ ___
    # ___ ___ ___
    # ___ ___ ___
    # ___ ___ ___
    #mean_init_z = make_init.makeinit(xv, yv, geoclaw_first) #My own function to generate the gaussian hump
    #If bowltopo not found, then use - 
    print "Creating the initial Gaussian hump ...\n"
    #mean_init_z = maketopo.qinit(xv, yv)
    mean_init_z = maketopo.qinit(xv/2.0, yv/2.0)

    #Create ensemble members based on the mean value vector
    #Ensemble data format
    # __ __ __ __ __ __ __ ... __
    # __ __ __ __ __ __ __ ... __
    # __ __ __ __ __ __ __ ... __
    # __ __ __ __ __ __ __ ... __
    # __ __ __ __ __ __ __ ... __
    # __ __ __ __ __ __ __ ... __
    # __ __ __ __ __ __ __ ... __
    # __ __ __ __ __ __ __ ... __
    print "Creating the first ensemble members ...\n"
    make_init_ens.makeinitens(mean_init_z,num_ens, "first_ens_")

    #for j in range(np.size(dtobs)-1):
    for k,j in enumerate(dtobs[:-1]):

        if DA:
            print "Data assimilation is on ...\n"
            print "##########################################################"
            print "----------------------Recevied observation at " + str(j) +"----------------------" 
            print "##########################################################"
        
        
            #Create observation data
            print "Creating observation data ...\n"
            if not firsttime:
                mean_init_z = np.loadtxt("ens_"+str((num_ens+1)/2)+".txt")
            observation = obs.make_obs(mx, my, dxobs, dyobs,stddev_obs,maketopo.qinit(mxv,myv), testing=False)
            
        #Write ensemble_tracker
        # For every forward run, Fortran reads the ensemble number from ens_tracker.
        with open("ens_tracker","w") as ens_tracker:
            ens_tracker.write("1")
    
        for i in range(1, num_ens+1):
            #Define pdaf input and output file names
            pdaf_input = "../ens_" + str(i) + ".txt"
            geoclaw_input = "hump_ens_" + str(i) + ".txt"
            pdaf_output = "../ens_0"+str(i)+"_ana.txt"
            #first_ensemble = "../first_ens_" + str(i) + ".txt"
            #read_geoclaw_output = "fort.q00" + str(output_times)
            subdir_name = "ens_"+str(i)

            print "\n##########################################################"
            print "------------- ensemble unit number " + str(i) + " at " + str(j) + " secs ----------" 
            print "##########################################################"
         
            #Check if path exists and create subfolders to run individual geoclaw
            cmdir.take(subdir_name+"_"+str(k))

            #Change to subdirectory
            os.chdir(subdir_name +"_"+str(k))

            #Prepare qinit files for geoclaw
            #Convert format of ensemble to data input format of Geoclaw qinit
            if firsttime:
                print "First time is executed\n" 
                #pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, first_ensemble, geoclaw_input)
                #pdaf_to_geoclaw.pdaf_to_geoclaw(xv/2.0, yv/2.0, first_ensemble, geoclaw_input)
                topotools.topo1writer(geoclaw_input,maketopo.qinit, xlower, xupper, ylower, yupper, nxpoints,nypoints)
            else:
                if DA:
                    pdaf_to_geoclaw.pdaf_to_geoclaw(mxv, myv, pdaf_output, geoclaw_input)
                else:
                    print "Second time is executed\n" 
                    #pdaf_to_geoclaw.pdaf_to_geoclaw(mxv, myv, pdaf_input, geoclaw_input)
                    np.savetxt("dummy1",np.zeros((mx,my)))
                    pdaf_to_geoclaw.pdaf_to_geoclaw(mxv, myv, "dummy1", geoclaw_input)
            
            #Run Geoclaw forecast step
            #---------------------------------------#
            ########        FORECAST       ##########
            #---------------------------------------#
            # Prepare individual geoclaw input
            print "Executing forecast step for ens_number" + str(k)
            run_geoclaw.run_geoclaw(k, dtobs=dtobs, mx = mx, my = my, geoclaw_input = geoclaw_input, topography = topo_path, output_times = output_times, max_amr = amr_max_level, geoclaw_exec = geoclaw_exec)
            #stop_case = ramr.ReadAmrForLevel(read_geoclaw_output, 1.0)

            
            #---------------------------------------------------------------#
            ########        GEOCLAW OUTPUT TO PDAF INPUT I/O       ##########
            #---------------------------------------------------------------#
            #np.savetxt(pdaf_input,reshaped_eta)
            #np.savetxt(pdaf_input,stop_case.eta_with_land)
            
            # Very dangerous
            if not (j == dtobs[-2]):
                shutil.copy2("fort.q0012","../fort.q0012_" + subdir_name)
                #pdb.set_trace()
            #Go back one directory  
            os.chdir("../")
        
        
        #Run PDAF assimilation step
        #-------------------------------------------#
        ########        ASSIMILATION       ##########
        #-------------------------------------------#
        if DA:
            subprocess.call([PDAF_executable])
        
        firsttime = False
        os.remove("ens_tracker")

    #-------------------------------------------#
    ###########     POST-PROCESSING    ##########
    #-------------------------------------------#
    #plot_ens.plot_ens(xv,yv,num_ens)
    if (num_ens != 1):
        plot_ens.plot_ens(xv,yv,num_ens, initial_state = True, analysis_state = False)
        
    #-----------------------------------------------------------------------
    #Verification of geoclaw restart
    #-----------------------------------------------------------------------
    verification_plot.verify_das(dtobs, (num_ens+1)/2, output_times, num_ens)

    if DA:
        plotmap.docontour(mxv,myv,observation, error_reshaped_masked_eta_land, "Observation data", -1000., 500)
        plotmap.docontour(mxv,myv,np.loadtxt("state_ana.txt"), error_reshaped_masked_eta_land,"Final state at last assimilation", -0.9, 0.9)


if __name__=='__main__':
    
    #Construct topography
    z1 = maketopo.maketopo()
    
    main()
