import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ReadAmrForLevel as ramr
import read_amr2
#import chunk
#import plotmap
#import myplots

def make_obs(mxv, myv, obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype):

    obs_x = np.linspace(xobs_start,xobs_end, nxobs).astype('int32')
    obs_y = np.linspace(yobs_start,yobs_end, nyobs).astype('int32')
    obs_xv, obs_yv = np.meshgrid(obs_x, obs_y)
    obs_mat = np.ones(np.shape(mxv))*-999.0

    for i,j in enumerate(obs_time_list):

        read_geoclaw_output = "./_output_original" + ictype + "/fort.q" + str(i+1).zfill(4)
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

        read_geoclaw_output = "../_output_original" + ictype + "/fort.q" + str(i+1).zfill(4)
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


def make_obs_testing2(obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype):

    obs_x = np.linspace(xobs_start,xobs_end, nxobs)
    obs_y = np.linspace(yobs_start,yobs_end, nyobs)
    #obs_x = np.array([1,5,7,15,19])
    #obs_y = np.array([3,7,15,19, 25])

    for i,j in enumerate(obs_time_list):

        read_geoclaw_output = "../_output_original" + ictype + "/fort.q" + str(i+1).zfill(4)
        original_case = read_amr2.ReadAmr2(read_geoclaw_output)
        original_case = original_case.amrdataframe()

        #orpd_df = orpd[orpd.xcoord.isin(obs_x) & orpd.ycoord.isin(obs_y)&(orpd.amrlevel == 1.0)]
        mask1 = np.isclose(original_case['xcoord'].values[:,None], obs_x).any(axis=1)
        mask2 = np.isclose(original_case['ycoord'].values[:,None], obs_y).any(axis=1)
        orpd_df = original_case[mask1 & mask2]
        #print orpd.xcoord
        #Construct observations
        obs_file = "obs_step"+str(j)+".txt"
        savefile = "obs_step" + str(j) + ".pdf"
        #print "Observation at chosen location - ",obs_mat[obs_xv,obs_yv]
        #print "Writing observation file - ", obs_file
        #print obs_mat
        #np.savetxt(obs_file, obs_mat, fmt = "%12.10f")
        orpd_df.to_csv(obs_file, sep = "\t", header=False, index=False, columns=["xcoord", "ycoord", "height"])

def make_obs_testing3(obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype):

    obs_x = np.linspace(xobs_start,xobs_end, nxobs)
    obs_y = np.linspace(yobs_start,yobs_end, nyobs)
    #obs_x = np.array([1,5,7,15,19])
    #obs_y = np.array([3,7,15,19, 25])

    for i,j in enumerate(obs_time_list):

        read_geoclaw_output = "../_output_main/" + ictype + "/fort.q" + str(j).zfill(4)
        print "Reading " + read_geoclaw_output
        original_case = read_amr2.ReadAmr2(read_geoclaw_output, sort_in_frame=False)
        original_case = original_case.pandas_dataframe

        #orpd_df = orpd[orpd.xcoord.isin(obs_x) & orpd.ycoord.isin(obs_y)&(orpd.amrlevel == 1.0)]
        mask1 = np.isclose(original_case['xcoord'].values[:,None], obs_x).any(axis=1)
        mask2 = np.isclose(original_case['ycoord'].values[:,None], obs_y).any(axis=1)
        orpd_df = original_case[mask1 & mask2]
        #print orpd.xcoord
        #Construct observations
        obs_file = "obs_step"+str(i+1)+".txt"
        savefile = "obs_step" + str(i+1) + ".pdf"
        print "Writing " + obs_file
        #print "Observation at chosen location - ",obs_mat[obs_xv,obs_yv]
        #print "Writing observation file - ", obs_file
        #print obs_mat
        #np.savetxt(obs_file, obs_mat, fmt = "%12.10f")
        orpd_df.to_csv(obs_file, sep = "\t", header=False, index=False, columns=["xcoord", "ycoord", "eta"])




if __name__=="__main__":
    #xobs_start = -119.94
    #yobs_start = -59.94
    #xobs_end = -71.94
    #yobs_end = -19.94
    xobs_start = -119.7
    yobs_start = -59.7
    xobs_end = -71.7
    yobs_end = -19.7
    nxobs = 21
    nyobs = 21
    ictype = ""
    #num_time_steps=6
    #obs_xrange = np.arange(10.0,11.0,1.0)
    #obs_yrange = np.arange(10.0,11.0,1.0)
    #obs_xmesh,obs_ymesh = np.meshgrid(obs_xrange, obs_yrange)
    #obs_time_list = np.linspace(180,780,11,dtype='int32')
    obs_time_list = np.linspace(10800/3600,46800/3600,11,dtype='int32')


    #make_obs(mxv, myv, np.linspace(20,200,10, dtype="int32"), xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype)
    ##make_obs_testing(mxv, myv, obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype)
    make_obs_testing3(obs_time_list, xobs_start, yobs_start, xobs_end, yobs_end, nxobs, nyobs, ictype)
    #if ((nx>50)&(ny>50)):
    #    chunk.chunk_write(obs_time_list,type1="obs")
