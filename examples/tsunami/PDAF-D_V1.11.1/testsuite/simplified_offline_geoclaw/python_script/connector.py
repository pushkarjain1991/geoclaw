import os
import sys
import subprocess
import shutil
import numpy as np
import os
import sys
import make_init_ens
import pdaf_to_geoclaw
import obs
import mesh_interpol
import cmdir
import ensemble_class 
import maketopo
import plot_ens
import scipy_interpol
import plotmap
from scipy import interpolate
import matplotlib.pyplot as plt
import read_amr as ramr


def main():
    """
    Performs forecast and assimilation steps of geoclaw
    """
    amr_max_level = 1
    output_times = 12
    DA = False
    num_ens = 9
    dtobs = [0.0,4.0, 8.0]
    
    
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
    
    #DA parameters
    stddev_obs = 0.5
    dxobs = 10
    dyobs = 10
    PDAF_executable = "./PDAF_offline"
   
    x = np.linspace(xlower,xupper,nxpoints)
    y = np.linspace(yupper,ylower,nypoints)
    xv,yv = np.meshgrid(x,y)
    
    dx = (xupper-xlower)/(nxpoints-1)
    dy = (yupper-ylower)/(nypoints-1)
    
    x_cell = np.linspace(xlower + dx/2.0, xupper - dx/2.0, mx)
    y_cell = np.linspace(ylower + dy/2.0, yupper - dy/2.0, my)
    mxv, myv = np.meshgrid(x_cell,y_cell)


    firsttime = True
    
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
            print "Data assimilation is on ...\n\n\n"
            print "##########################################################"
            print "----------------------Recevied observation at " + str(j) +"----------------------" 
            print "##########################################################"
        
        
            #Create observation data
            print "Creating observation data ...\n"
            observation = obs.make_obs(nxpoints, nypoints, dxobs, dyobs,stddev_obs,mean_init_z, testing=False)
    
        for i in range(1, num_ens+1):
            #Define pdaf input and output file names
            pdaf_input = "../ens_" + str(i) + ".txt"
            geoclaw_input = "hump_ens_" + str(i) + ".txt"
            pdaf_output = "../ens_0"+str(i)+"_ana.txt"
            first_ensemble = "../first_ens_" + str(i) + ".txt"
            topo_path = "../bowl.topotype2"

            print "##########################################################"
            print "------------- ensemble unit number " + str(i) + " at " + str(j) + " secs ----------" 
            print "##########################################################"
         
            #Check if path exists and create subfolders to run individual geoclaw
            subdir_name = "ens_"+str(i)
            cmdir.take(subdir_name)

            #Change to subdirectory
            os.chdir(subdir_name)

            #Prepare qinit files for geoclaw
            #Convert format of ensemble to data input format of Geoclaw qinit
            if firsttime:
                #pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_input, geoclaw_input)
                #This is an unnecessary step. Individual z vectors
                # are already calculated by makeinitens. 
                print "\n\n\n\n\nFirst time is executed\n\n\n\n\n" 
                #pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, first_ensemble, geoclaw_input)
                pdaf_to_geoclaw.pdaf_to_geoclaw(xv/2.0, yv/2.0, first_ensemble, geoclaw_input)
            else:
                if DA:
                    pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_output, geoclaw_input)
                else:
                    print "\n\n\n\n\nSecond time is executed\n\n\n\n\n" 
                    pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_input, geoclaw_input)
            
            #Run Geoclaw forecast step
            #---------------------------------------#
            ########        FORECAST       ##########
            #---------------------------------------#
            # Prepare individual geoclaw input
            #subprocess.call(["make",".output"])
            hello = ensemble_class.ensemble()
            hello.rundata.clawdata.t0 = dtobs[k]
            hello.rundata.clawdata.tfinal = dtobs[k+1]
            hello.rundata.qinit_data.qinitfiles[-1][-1]=geoclaw_input
            #hello.rundata.qinit_data.qinitfiles[-1]=[1,2,geoclaw_input]
            #hello.rundata.topo_data.topofiles[-1]=[2, 1, 1, 0., 1.e10, topo_path]
            hello.rundata.topo_data.topofiles[-1][-1]=topo_path
            hello.rundata.clawdata.num_output_times = output_times
            hello.rundata.amrdata.amr_levels_max=amr_max_level
            hello.rundata.write()
            print "yoyoyoyo\n\n\n"
            subprocess.call(geoclaw_exec)

            #---------------------------------------------------------------#
            ########        GEOCLAW OUTPUT TO PDAF INPUT I/O       ##########
            #---------------------------------------------------------------#
            #Extract water surface elevation from geoclaw fort.q file
            read_geoclaw_output = "fort.q00" + str(output_times)

            amrread = ramr.ReadAmr(read_geoclaw_output)
            amrframe = amrread.amrdataframe()

            #total_height, eta = np.loadtxt(read_geoclaw_output, skiprows=9, usecols = (0,3), unpack=True)
            total_height1 = amrframe["height"][amrframe.amrlevel == 1.0] 
            eta1 = amrframe["eta"][amrframe.amrlevel == 1.0] 
            if (i == (num_ens+1)/2):
                #np.savetxt("../eta_with_land_" + str(j) + ".txt", eta)
                eta1.to_csv("../eta_with_land_" + str(j) + ".txt", sep=' ', index=False)
            eta1[total_height1 == 0.0] = 0.0
            if (i == (num_ens+1)/2):
                #np.savetxt("../eta_before_interp_" + str(j) + ".txt", eta)
                eta1.to_csv("../eta_before_interp_" + str(j) + ".txt", sep=' ', index=False)
        
            eta = eta1.as_matrix()
            #Write new ensembles into PDAF input format
            reshaped_eta = np.reshape(eta,(mx,my))
            #sys.exit(0)
            #Interpolate eta values from cell centers to nodes. interp_eta will be of size nx*ny
            #Will be a problem if AMR starts. 
            #Option 1
            interp_eta = mesh_interpol.interpol(reshaped_eta)

            #Option 2
            #f = interpolate.interp2d(mxv, myv, reshaped_eta)
            #interp_eta = f(x,y)

            #Option 3
           #interp_eta = scipy_interpol.interpol(mxv,myv,reshaped_eta)

            if (i == (num_ens+1)/2):
                np.savetxt("../eta_after_interp_" + str(j) + ".txt", np.reshape(interp_eta, (nxpoints*nypoints,1)))
            
            #print np.shape(interp_eta)
            np.savetxt(pdaf_input,interp_eta)
            
            # Very dangerous
            if not (j == dtobs[-2]):
                print j
                shutil.copy2("fort.q0012","../")

            #Go back one directory  
            os.chdir("../")
        
        
        #Run PDAF assimilation step
        #-------------------------------------------#
        ########        ASSIMILATION       ##########
        #-------------------------------------------#
        if DA:
            subprocess.call([PDAF_executable])
        
        firsttime = False

        #-------------------------------------------#
        ###########     POST-PROCESSING    ##########
        #-------------------------------------------#
    #plotmap.docontour(x_cell,y_cell,reshaped_eta)
    #plotmap.docontour(xv,yv,interp_eta)
    #plotmap.docontour(x,y,interp_eta)
    #plot_ens.plot_ens(xv,yv,num_ens)
    #plot_ens.plot_ens(xv,yv,num_ens, initial_state = True, analysis_state = False)
    plotmap.docontour(mxv,myv,reshaped_eta)
    plotmap.docontour(xv,yv,interp_eta)
    if DA:
        plotmap.docontour(xv,yv,observation)
        plotmap.docontour(xv,yv,np.loadtxt("state_ana.txt"))
    


#Just end everything and move on ahead in life


if __name__=='__main__':
    
    #Construct topography
    z1 = maketopo.maketopo()
    
    main()
