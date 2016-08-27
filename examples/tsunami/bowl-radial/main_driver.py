import numpy as np
from inputs_online import make_obs
from inputs_online import error_analysis
from inputs_online import create_mask
from inputs_online import run_geoclaw
from inputs_online import make_init_ens
from inputs_online import error_analysis
from inputs_online import chunk

if __name__ == "__main__":
    
    #------------------------------------------------------#
    #Follow these instructions to run GeoCLAW manually
    #1. Use -cpp as only #FFLAGS in Makefile
    #2. Run "make .output"
    #3. Move the output to _output_original_hump. This output will serve as observation values for assimilation
    #------------------------------------------------------#
    
    nx = 100
    ny = 100
    xlower = -99.0
    xupper = 99.0
    ylower = -99.0
    yupper = 99.0
    num_ens = 1
    ictype = "hump"
    total_steps = 200
    output_step_interval = 20
    delt_obs = output_step_interval
    #rms_obs = 0.0000001
    rms_obs = 0.00000001
    filtertype = 2
    
    obs_time_list = np.linspace(output_step_interval, total_steps,total_steps/output_step_interval, dtype='int32')
    num_time_steps = total_steps/output_step_interval + 1
    X = np.linspace(xlower, xupper, nx)
    Y = np.linspace(ylower, yupper, ny)
    xv,yv = np.meshgrid(X,Y)
   
    #Make initial ensemble
    make_init_ens.makeinitens2(xv, yv, num_ens,ictype = ictype)
    #make_init_ens.makeinitens(xv, yv, num_ens,ictype=ictype)
    if ((nx > 50) & (ny > 50)):
        chunk.chunk_write(range(1,num_ens+1),type1="init")
   
    #Make observations
    #if num_ens != 1:
    xobs_start = 20 
    yobs_start = 20
    xobs_end = 80
    yobs_end = 80
    nxobs = 4
    nyobs = 4
    make_obs.make_obs(xv, yv, obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype)
        #make_obs.make_obs_testing(xv, yv, obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype)
    if ((nx > 50) & (ny > 50)):
        chunk.chunk_write(obs_time_list,type1="obs")

    #------------------------------------------------------#
    #Follow these instructions to run assimilation version of GeoCLAW
    #1. Use the alternate #FFLAGS in Makefile
    #2. Make sure setrun.py has num_ens set to same value as above num_ens
    #3. Recompile GeoCLAW after "make clobber"
    #4. Run "make .output"
    #5. Proceed to error analysis
    #------------------------------------------------------#
    run_geoclaw.run_geoclaw("../xgeoclaw",\
            filtertype=filtertype, \
            num_ens=num_ens, \
            rms_obs=rms_obs,\
            delt_obs=delt_obs,\
            output_step_interval=output_step_interval,\
            total_steps=total_steps)

    #------------------------------------------------------#
    #Error analysis
    error_analysis.error_calc(ictype, num_time_steps, ploterror=True)
    #------------------------------------------------------#

    #------------------------------------------------------#
    #Vizualization
    #------------------------------------------------------#
    #Create mask of land and water
    #zmin = 80.
    #create_mask(xv,yv,zmin)

