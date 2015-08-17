import os
import sys
import subprocess
import shutil
import numpy as np
import os
import sys
import make_init_ens
import pdaf_to_geoclaw
import remove_file
import obs
import mesh_interpol
import cmdir
import ensemble_class 
import maketopo
import plot_ens2
import scipy_interpol
import plotmap
from scipy import interpolate
import matplotlib.pyplot as plt


def main():
    """
    Performs forecast and assimilation steps of geoclaw
    """
    #Model Parameters
    nxpoints = 51
    mx = nxpoints-1
    nypoints = 51
    my = nypoints-1
    #xlower = -50.e0
    #xupper = 50.e0
    #yupper = 50.e0
    #ylower = -50.e0
    xlower = -100.e0
    xupper = 100.e0
    yupper = 100.e0
    ylower = -100.e0
    geoclaw_exec = "../xgeoclaw"
    amr = False
    
    #DA parameters
    DA = False
    num_ens = 1
    stddev_obs = 0.5
    dxobs = 5
    dyobs = 4
    #dtobs = [2.0, 4.0, 6.0, 8.0]
    #dtobs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    #dtobs = [0.0, 2.0, 4.0, 6.0, 8.0]
    dtobs = [0.0, 4.0, 8.0]
    #dtobs = [0.0, 4.0]
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
    mean_init_z = maketopo.qinit(xv, yv)

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
        print k,j
        print "yoyoyoy3\n\n\n"

        if DA == True:
            print "Data assimilation is on ...\n\n\n"
            print "##########################################################"
            print "----------------------Recevied observation at " + str(j) +"----------------------" 
            print "##########################################################"
        
        
            #Create observation data
            print "Creating observation data ...\n"
            obs.make_obs(nxpoints, nypoints, dxobs, dyobs,stddev_obs,mean_init_z)
    
        for i in range(1, num_ens+1):
            print "yoyoyoy5\n\n\n"
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
            if firsttime == True:
                #print "yoyoyoy6\n\n\n"
                #pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_input, geoclaw_input)
                #This is an unnecessary step. Individual z vectors
                # are already calculated by makeinitens. 
                print "\n\n\n\n\nFirst time is executed\n\n\n\n\n" 
                pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, first_ensemble, geoclaw_input)
            else:
                if DA == True:
                    #print "yoyoyoy7\n\n\n"
                    pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_output, geoclaw_input)
                else:
                    #print "yoyoyoy8\n\n\n"
                    print "\n\n\n\n\nSecond time is executed\n\n\n\n\n" 
                    pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_input, geoclaw_input)
            
            #print "yoyoyoy9\n\n\n"
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
            hello.rundata.clawdata.num_output_times = 12
            #hello.rundata.clawdata.num_output_times = 24
            if amr == False:
                hello.rundata.amrdata.amr_levels_max=1
            else:
                hello.rundata.amrdata.amr_levels_max=2
            hello.rundata.write()
            subprocess.call(geoclaw_exec)

            #Extract water surface elevation from geoclaw fort.q file
                #Extract land nodes
            read_geoclaw_output = "fort.q00" + str(hello.rundata.clawdata.num_output_times)
            total_height, eta = np.loadtxt(read_geoclaw_output, skiprows=9, usecols = (0,3), unpack=True)
            if (i == (num_ens+1)/2):
                np.savetxt("../eta_with_land_" + str(j) + ".txt", eta)
            getland = [total_height == 0.0]
            eta[getland]= 0.0
            #np.savetxt("myeta" + str(j),eta)
            if (i == (num_ens+1)/2):
                np.savetxt("../eta_before_interp_" + str(j) + ".txt", eta)
                #z1 = np.loadtxt("../yoyohoney", unpack=True, usecols=[0])
                #print z1
                #np.reshape(z1, (50,50))
                #print z1
        
            #Write new ensembles into PDAF input format
            reshaped_eta = np.reshape(eta,(mx,my))
            print reshaped_eta
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
        
            #Go back one directory  
            os.chdir("../")
        
        
        #Run PDAF assimilation step
        #-------------------------------------------#
        ########        ASSIMILATION       ##########
        #-------------------------------------------#
        print "yoyoyoy10\n\n\n"
        if DA == True:
            print "yoyoyoy11\n\n\n"
            subprocess.call([PDAF_executable])
        
        firsttime = False
        print "yoyoyoy12\n\n\n"

        #-------------------------------------------#
        ###########     POST-PROCESSING    ##########
        #-------------------------------------------#
    #plotmap.docontour(x_cell,y_cell,reshaped_eta)
    #plotmap.docontour(xv,yv,interp_eta)
        #plot_ens2.plot_ens(xv,yv,num_ens)
        #plot_ens2.plot_ens(xv,yv,num_ens, initial_state = "True", analysis_state = "False")
    #plotmap.docontour(x_cell,y_cell,reshaped_eta)
    plotmap.docontour(np.linspace(-100.0 + 1.0, 100.0 - 1.0, mx),np.linspace(-100.0 + 1.0, 100.0 - 1.0, my),reshaped_eta)
    plotmap.docontour(np.linspace(-100.0, 100.0, nxpoints),np.linspace(-100.0,100.0,nypoints),interp_eta)
    


#Just end everything and move on ahead in life


if __name__=='__main__':
    
    #Construct topography
    z1 = maketopo.maketopo()
    
    main()
