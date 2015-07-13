import os
import sys
import subprocess
import shutil
import numpy as np
import os
import sys
import make_init
import make_init_ens
import pdaf_to_geoclaw
import remove_file
import obs
import mesh_interpol
import cmdir
import ensemble_class 
import geoclaw_input_format as gcif


# First will try running using qinit.

#Create initial ensemble

#----------------------------------------------------#
#Do the forecast step
#----------------------------------------------------#
#Use the data of ens_n.txt to re-initialize the data.
#Note that the output is written to fort.q. The 4th column is the eta
# Write it in hump.xyz format [x y val]

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
    
    #DA parameters
    num_ens = 9
    obs_t_interval = 10
    stddev_obs = 0.5
    dxobs = 5
    dyobs = 4
    #dtobs = [2.0, 4.0, 6.0, 8.0]
    dtobs = [2.0, 4.0, 6.0, 8.0]
    firsttime = True
    PDAF_executable = "./PDAF_offline"
    geoclaw_exec = "./xgeoclaw"

    geoclaw_input = "hump.xyz"     
    radialbowl_files = ["bowl.topotype2", "Makefile", "setrun.py", geoclaw_exec]
    radialbowl_path = "/h2/pkjain/Desktop/Pushkar/clawpack/geoclaw/examples/tsunami/bowl-radial/"
   
    x = np.linspace(xlower,xupper,nxpoints)
    y = np.linspace(yupper,ylower,nypoints)
    xv,yv = np.meshgrid(x,y)
        
    #Create main initial data file in format of the Geoclaw input
    #Create hump.xyz
    mean_init_z = make_init.makeinit(xv, yv, geoclaw_input)

    #Create ensemble members based on the mean value vector
    make_init_ens.makeinitens(mean_init_z,num_ens, "ens_")
    
    for j in range(np.size(dtobs)-1):
    #for j in range(2):
        #Create observation data
        obs.make_obs(nxpoints, nypoints, dxobs, dyobs,stddev_obs,mean_init_z)

    
        #Convert format of ensemble to data input format of Geoclaw qinit
        for i in range(1, num_ens+1):

            print "##########################################################"
            print "----------------------" + str(i) + "----------------------" 
            print "##########################################################"
         
            #Check if path exists and create subfolders to run individual geoclaw
            subdir_name = "ens_"+str(i)
            cmdir.take(subdir_name)

            #Change to subdirectory
            os.chdir(subdir_name)

            #Define pdaf input and output file names
            pdaf_input = "../ens_" + str(i) + ".txt"
            pdaf_output = "../ens_0"+str(i)+"_ana.txt"

            #Prepare qinit files for geoclaw
            if firsttime:
                #pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_input, geoclaw_input)
                #This is an unnecessary step. Indiavidual z vetors
                # are already calculated by makeinitens. 
                z_ind = np.loadtxt(pdaf_input)
                print z_ind
                gcif.geoclaw_input_format(xv, yv, z_ind,pdaf_input)
                #firsttime = False
            else:
                pdaf_to_geoclaw.pdaf_to_geoclaw(xv, yv, pdaf_output, geoclaw_input)

            #Copy the radial bowl test case files to every sub directory
            for files in radialbowl_files:
                shutil.copy2(os.path.join(radialbowl_path,files),os.getcwd())
        
            #---------------------------------------#
            ########        FORECAST       ##########
            #---------------------------------------#
            #Run Geoclaw forecast step
            #subprocess.call(["make",".output"])
            hello = ensemble_class.ensemble()
            hello.rundata.clawdata.t0 = dtobs[j]
            hello.rundata.clawdata.tfinal = dtobs[j+1]
            hello.rundata.qinit_data.qinitfiles[-1]=[1,2,"../ens_"+str(i)+".txt"]
            #print hello.rundata.qinit_data.qinitfiles
            hello.rundata.write()
            subprocess.call(geoclaw_exec)

            #Extract water surface elevation from geoclaw fort.q file
            eta = np.loadtxt("fort.q00" + str(obs_t_interval), skiprows=9, usecols = [3])
            #np.savetxt("../ens_"+str(i)+".txt_new", eta)

        
            #Write new ensembles into PDAF input format
            reshaped_eta = np.reshape(eta,(mx,my))
            #print np.shape(eta)
            #Interpolate eta values from cell centers to nodes. interp_eta will be of size nx*ny
            interp_eta = mesh_interpol.interpol(reshaped_eta)
            #print np.shape(interp_eta)
        
            #Go back one directory  
            os.chdir("../")
        
            np.savetxt("../ens_"+str(i)+".txt",interp_eta)
            #np.savetxt("../ens_"+str(i)+".txt_new_reshaped",interp_eta)
        
        #Run PDAF assimilation step
        #-------------------------------------------#
        ########        ASSIMILATION       ##########
        #-------------------------------------------#
        subprocess.call([PDAF_executable])
        
        # Run Geoclaw forecast step for all the ensemble members
        
        firsttime=False


 

#Do the assimilation step


#Just end everything and move on ahead in life


if __name__=='__main__':
    main()
