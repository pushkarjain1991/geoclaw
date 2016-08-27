import plotmap
import ReadAmrForLevel as ramr
import myplots
import numpy as np
import matplotlib.pyplot as plt
import read_amr
import chunk


def make_obs(mxv, myv, obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype):

    obs_x = np.linspace(xobs_start,xobs_end, nxobs).astype('int32')
    obs_y = np.linspace(yobs_start,yobs_end, nyobs).astype('int32')
    obs_xv, obs_yv = np.meshgrid(obs_x, obs_y)
    obs_mat = np.ones(np.shape(mxv))*-999.0

    for i,j in enumerate(obs_time_list):

        read_geoclaw_output = "./_output_original_" + ictype + "/fort.q" + str(i+1).zfill(4)
        original_case = ramr.ReadAmrForLevel(read_geoclaw_output, 1.0)
        original_water = original_case.water
        original_land = original_case.land
        #mxv = original_case.mxv
        #myv = original_case.myv

        #Construct observations
        obs_file = "obs_step"+str(j)+".txt"
        savefile = "obs_step" + str(j) + ".pdf"
        obs_mat[obs_xv,obs_yv] = original_water[obs_xv,obs_yv]

        print "Observation at chosen location - ",obs_mat[obs_xv,obs_yv] 
        print "Writing observation file - ", obs_file
        np.savetxt(obs_file, obs_mat, fmt = "%12.10f")
        obs_mat_water = np.ma.array(obs_mat,mask = original_case.land==0.0)
        obs_mat_land = original_case.land


        plotmap.docontour(mxv,myv,obs_mat_water, obs_mat_land, "Observation", -999.0, 1.0, savefile=savefile)
    

def make_obs_testing(mxv, myv, obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype):

    obs_x = np.linspace(xobs_start,xobs_end, nxobs).astype('int32')
    obs_y = np.linspace(yobs_start,yobs_end, nyobs).astype('int32')
    obs_xv, obs_yv = np.meshgrid(obs_x, obs_y)
    obs_mat = np.ones((100,100))*-999.0

    for i,j in enumerate(obs_time_list):

        read_geoclaw_output = "../_output_original_" + ictype + "/fort.q" + str(i+1).zfill(4)
        #original_case = ramr.ReadAmrForLevel(read_geoclaw_output, 1.0)
        #original_water = original_case.water
        #original_land = original_case.land
        #mxv = original_case.mxv
        #myv = original_case.myv

        original_case = read_amr.ReadAmr(read_geoclaw_output)
        orpd = original_case.pandas_dataframe
        #print orpd["eta"][(orpd.xcoord == -99.0)&(orpd.ycoord == 99.0)&(orpd.amrlevel == 1.0)]
        for i1 in obs_x:
            for j1 in obs_y:
                obs_mat[i1,j1] = orpd["eta"][(orpd.xcoord == mxv[i1,j1])&(orpd.ycoord == myv[i1,j1])&(orpd.amrlevel == 1.0)]
        #print orpd["eta"][(orpd.xcoord == np.ravel(mxv[obs_xv]))&(orpd.ycoord == np.ravel(myv[obs_yv]))&(orpd.amrlevel == 1.0)]
        #Construct observations
        obs_file = "obs_step"+str(j)+".txt"
        savefile = "obs_step" + str(j) + ".pdf"
        print "Observation at chosen location - ",obs_mat[obs_xv,obs_yv] 
        print "Writing observation file - ", obs_file
        np.savetxt(obs_file, obs_mat, fmt = "%12.10f")
    

if __name__=="__main__":
    nx = 100
    ny = 100
    x = np.linspace(-99,99,nx)
    y = np.linspace(-99,99,ny)
    mxv,myv = np.meshgrid(x,y)
    #print mxv[24][30]
    xobs_start = 45 
    yobs_start = 45
    xobs_end = 55
    yobs_end = 55
    nxobs = 5
    nyobs = 5
    ictype = "hump"
    #num_time_steps=6
    #obs_xrange = np.arange(10.0,11.0,1.0)
    #obs_yrange = np.arange(10.0,11.0,1.0)
    #obs_xmesh,obs_ymesh = np.meshgrid(obs_xrange, obs_yrange)
    obs_time_list = np.linspace(20,200,10,dtype='int32')


    make_obs(mxv, myv, np.linspace(20,200,10, dtype="int32"), xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype)
    #make_obs_testing(mxv, myv, obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype)
    if ((nx>50)&(ny>50)):
        chunk.chunk_write(obs_time_list,type1="obs")
