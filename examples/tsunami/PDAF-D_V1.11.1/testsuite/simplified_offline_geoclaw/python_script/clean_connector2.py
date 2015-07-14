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


def main():
    """
    Performs forecast and assimilation steps of geoclaw
    """
    #Model Parameters
    nxpoints = 51
    mx = nxpoints-1
    nypoints = 51
    my = nypoints-1
    xlower = -50.e0
    xupper = 50.e0
    yupper = 50.e0
    ylower = -50.e0
    geoclaw_exec = "../xgeoclaw"
    
    #DA parameters
    num_ens = 9
    stddev_obs = 0.5
    dxobs = 5
    dyobs = 4
    dtobs = [2.0, 4.0, 6.0, 8.0]
    firsttime = True
    
    PDAF_executable = "./PDAF_offline"
   
    x = np.linspace(xlower,xupper,nxpoints)
    y = np.linspace(yupper,ylower,nypoints)
    xv,yv = np.meshgrid(x,y)
    
    #Create main initial data file in format of the Geoclaw input
    #Create hump.xyz
    #Geoclaw input format
    # ___ ___ ___
    # ___ ___ ___
    # ___ ___ ___
    # ___ ___ ___
    # ___ ___ ___
    # ___ ___ ___
    #mean_init_z = make_init.makeinit(xv, yv, geoclaw_first)
    #If bowltopo not found, then use - 
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
    make_init_ens.makeinitens(mean_init_z,num_ens, "first_ens_")
    
    for j in range(np.size(dtobs)-1):
        print "##########################################################"
        print "----------------------Recevied observation at " + str(dtobs[j]) +"----------------------" 
        print "##########################################################"
        
        
        #Create observation data
        obs.make_obs(nxpoints, nypoints, dxobs, dyobs,stddev_obs,mean_init_z)

    
        for i in range(1, num_ens+1):
            #Define pdaf input and output file names
            pdaf_input = "../ens_" + str(i) + ".txt"
            geoclaw_input = "hump_+ens_" + str(i) + ".txt"
            pdaf_output = "../ens_0"+str(i)+"_ana.txt"
            first_ensemble = "../first_ens_" + str(i) + ".txt"
            topo_path = "../bowl.topotype2"

            print "##########################################################"
            print "---------------------- ensemble unit number " + str(i) + "----------------------" 
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
                pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, first_ensemble, geoclaw_input)
            else:
                pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_output, geoclaw_input)

            #---------------------------------------#
            ########        FORECAST       ##########
            #---------------------------------------#
            #Run Geoclaw forecast step
            #subprocess.call(["make",".output"])
            hello = ensemble_class.ensemble()
            hello.rundata.clawdata.t0 = dtobs[j]
            hello.rundata.clawdata.tfinal = dtobs[j+1]
            hello.rundata.qinit_data.qinitfiles[-1]=[1,2,geoclaw_input]
            hello.rundata.topo_data.topofiles[-1]=[2, 1, 1, 0., 1.e10, topo_path]
            hello.rundata.write()
            subprocess.call(geoclaw_exec)

            #Extract water surface elevation from geoclaw fort.q file
            eta = np.loadtxt("fort.q00" + str(hello.rundata.clawdata.num_output_times), skiprows=9, usecols = [3])

        
            #Write new ensembles into PDAF input format
            reshaped_eta = np.reshape(eta,(mx,my))
            #Interpolate eta values from cell centers to nodes. interp_eta will be of size nx*ny
            interp_eta = mesh_interpol.interpol(reshaped_eta)
            
            np.savetxt(pdaf_input,interp_eta)
        
            #Go back one directory  
            os.chdir("../")
        
        
        #Run PDAF assimilation step
        #-------------------------------------------#
        ########        ASSIMILATION       ##########
        #-------------------------------------------#
        subprocess.call([PDAF_executable])
        
        firsttime=False


#Just end everything and move on ahead in life


if __name__=='__main__':
    
    #Construct topography
    z1 = maketopo.maketopo()
    
    main()
