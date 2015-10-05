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
import read_amr
import ReadAmrForLevel as ramr


def main():
    """
    Performs forecast and assimilation steps of geoclaw
    """
    amr_max_level = 1
    output_times = 12
    DA = False
    num_ens = 1
    dtobs = [0.0, 2.0, 4.0, 6.0, 8.0]
    #dtobs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
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
            if not firsttime:
                mean_init_z = np.loadtxt("ens_"+str((num_ens+1)/2)+".txt")
            #observation = obs.make_obs(nxpoints, nypoints, dxobs, dyobs,stddev_obs,mean_init_z, testing=False)
            observation = obs.make_obs(mx, my, dxobs, dyobs,stddev_obs,mean_init_z, testing=False)
            
        #Write ensemble_tracker
        # For every forward run, Fortran reads the ensemble number from ens_tracker.
        with open("ens_tracker","w") as ens_tracker:
            ens_tracker.write("1")
    
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
                    pdaf_to_geoclaw.pdaf_to_geoclaw(mxv, myv, pdaf_output, geoclaw_input)
                else:
                    print "\n\n\n\n\nSecond time is executed\n\n\n\n\n" 
                    pdaf_to_geoclaw.pdaf_to_geoclaw(mxv, myv, pdaf_input, geoclaw_input)
            
            #Run Geoclaw forecast step
            #---------------------------------------#
            ########        FORECAST       ##########
            #---------------------------------------#
            # Prepare individual geoclaw input
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
            subprocess.call(geoclaw_exec)

            #---------------------------------------------------------------#
            ########        GEOCLAW OUTPUT TO PDAF INPUT I/O       ##########
            #---------------------------------------------------------------#
            #Extract water surface elevation from geoclaw fort.q file
            read_geoclaw_output = "fort.q00" + str(output_times)

            #amrread = ramr.ReadAmr(read_geoclaw_output)
            ##amrframe = amrread.amrdataframe()

            #total_height1 = amrread.get_mycolumn("height",amrl=1.0)
            #momx = amrread.get_mycolumn("xvel",amrl=1.0)
            #eta1 = amrread.get_mycolumn("eta",amrl=1.0)
            #eta1[total_height1 == 0.0] = 0.0
            #eta = eta1.as_matrix()
            #reshaped_eta = np.reshape(eta,(mx,my))
            
            ##Land
            #masked_eta_land = np.ma.array(total_height1, mask=momx==0.0E0)
            #reshaped_masked_eta_land = np.reshape(masked_eta_land,(mx,my))
            ##np.savetxt("../land",reshaped_masked_eta_land)
            
            ##Water
            #masked_eta_water = np.ma.array(eta1,mask=total_height1==0.0E0)
            #reshaped_masked_eta_water = np.reshape(masked_eta_water,(mx,my))
            ##np.savetxt("../water",reshaped_masked_eta_water)

            stop_case = ramr.ReadAmrForLevel(read_geoclaw_output, 1.0)

            #-----------------------------------------------------------------------
            #Verify with original data
            #-----------------------------------------------------------------------
            #original_fortq_file = str(4*output_times*int(dtobs[k+1])/int(dtobs[-1])).zfill(2)
            original_fortq_file = str(2*output_times*int(dtobs[k+1])/int(dtobs[-1])).zfill(2)
            print "Original file read is fort.q00" + original_fortq_file
            
            #original_case = ramr.ReadAmrForLevel("../original_radial_bowl/fort.q00" + original_fortq_file, 1.0)
            original_case = ramr.ReadAmrForLevel("../original_radial_bowl_48/fort.q00" + original_fortq_file, 1.0)
            
            #original_read = ramr.ReadAmr("../original_radial_bowl/fort.q00" + original_fortq_file)
            ##original_read = ramr.ReadAmr("../original_radial_bowl_48/fort.q00" + original_fortq_file)
            ##original_amrframe = original_read.amrdataframe()
            #original_total_height1 = original_read.get_mycolumn("height", amrl=1.0) 
            #original_momx = original_read.get_mycolumn("xvel", amrl = 1.0)
            #original_eta1 = original_read.get_mycolumn("eta",amrl=1.0)
            #original_eta1[original_total_height1 == 0.0] = 0.0
            #original_eta = original_eta1.as_matrix()
            #reshaped_original_eta = np.reshape(original_eta,(mx,my))

            #error_eta_percent = (eta - original_eta)*100.0/original_eta
            error_eta_percent = (stop_case.eta_with_land - original_case.eta_with_land)*100.0/original_case.eta_with_land
            error_reshaped_eta_percent = np.reshape(error_eta_percent,(mx,my))

            #error_eta = (eta - original_eta) 
            error_eta = (stop_case.eta_with_land - original_case.eta_with_land)
            error_reshaped_eta = np.reshape(error_eta,(mx,my))

            #Land _Original
            #error_masked_eta_land = np.ma.array(original_total_height1, mask=original_momx==0.0E0)
            #error_reshaped_masked_eta_land = np.reshape(error_masked_eta_land,(mx,my))
            error_reshaped_masked_eta_land = original_case.land
            #np.savetxt("../land",reshaped_masked_eta_land)
            
            #Water _Original
            error_masked_eta_water = np.ma.array(error_eta,mask=original_case.total_height==0.0E0)
            error_masked_eta_water_percent = np.ma.array(error_eta_percent,mask=original_case.total_height==0.0E0)
            error_reshaped_masked_eta_water = np.reshape(error_masked_eta_water,(mx,my))
            error_reshaped_masked_eta_water_percent = np.reshape(error_masked_eta_water_percent,(mx,my))
            #np.savetxt("../water",reshaped_masked_eta_water)

            if (i == (num_ens+1)/2):
                 plotmap.docontour(mxv,myv,stop_case.water, stop_case.land, "WSE: ens_number = "+str(i)+"; num_ens = " + str(num_ens)+"; time = " + str(dtobs[k+1]), -0.9, 0.9)
                 plotmap.docontour(mxv,myv,error_reshaped_masked_eta_water, error_reshaped_masked_eta_land, "WSE error: ens_number = "+str(i)+"; num_ens = " + str(num_ens)+"; time = " + str(dtobs[k+1]), -0.9,0.9)
                 plotmap.docontour(mxv,myv,error_reshaped_masked_eta_water_percent, error_reshaped_masked_eta_land, "WSE error %: ens_number = "+str(i)+"; num_ens = " + str(num_ens)+"; time = " + str(dtobs[k+1]), -200.0,200.0)
                 #plotmap.docontour(mxv,myv,error_reshaped_masked_eta_water, error_reshaped_masked_eta_land, "WSE error: ens_number = "+str(i)+"; num_ens = " + str(num_ens)+"; time = " + str(dtobs[k+1]), -0.25,0.17)
            
            #np.savetxt(pdaf_input,reshaped_eta)
            np.savetxt(pdaf_input,stop_case.eta_with_land)
            
            # Very dangerous
            if not (j == dtobs[-2]):
                #print j
                shutil.copy2("fort.q0012","../fort.q0012_" + subdir_name)

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
    #plot_ens.plot_ens(xv,yv,num_ens)
    if (num_ens != 1):
        plot_ens.plot_ens(xv,yv,num_ens, initial_state = True, analysis_state = False)
        
    #plotmap.docontour(xv,yv,interp_eta, "Interpolated WSE")
    if DA:
        plotmap.docontour(xv,yv,observation, "Observation data")
        plotmap.docontour(xv,yv,np.loadtxt("state_ana.txt"), "Final state at last assimilation")
    


#Just end everything and move on ahead in life


if __name__=='__main__':
    
    #Construct topography
    z1 = maketopo.maketopo()
    
    main()
