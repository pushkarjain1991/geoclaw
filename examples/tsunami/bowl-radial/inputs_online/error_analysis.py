import numpy as np
import pandas as pd
from scipy.linalg import norm
import matplotlib.pyplot as plt
import plotmap
import ReadAmrForLevel as ramr
import finderror


def calcerror(ictype, num_time_steps):
    for t in range(num_time_steps):
        original_file = "../_output_original_" + ictype + "/fort.q000" + str(t)
        df_original = pd.read_csv(original_file, dtype = np.float64, header=None, names=["total_height", "xmom", "ymom", "eta"], index_col=False, skip_blank_lines=True, skiprows=8, sep=r"\s+")

        test_file = "../_output/fort.q000" + str(t)
        df_test = pd.read_csv(test_file, dtype = np.float64, header=None, names=["total_height", "xmom", "ymom", "eta"], index_col=False, skip_blank_lines=True, skiprows=8, sep=r"\s+")

        #error = df_test["eta"] - df_original["eta"]
        error = df_test["total_height"] - df_original["total_height"]
        #error_pct = error/df_original["eta"]
        error_pct = error*100.0/df_original["total_height"]
        
        print "timestep = " + str(t)
        print "L_2 norm = ", norm(error) 
        print "L_inf norm = ", norm(error,np.inf)
        print "Max error percentage = ",norm(error_pct[df_original["total_height"]!=0.0],np.inf)
    
    return error, error_pct

def print_full(x, filename):
    pd.set_option('display.max_rows', len(x))
    #print(x)
    pd.reset_option('display.max_rows')
    x.to_csv(filename,sep='\t')

def error_calc(ictype, num_time_steps, ploterror=True):
    for i in range(num_time_steps):
        original_file = "../_output_original_" + ictype + "/fort.q" + str(i).zfill(4)
        original_case = ramr.ReadAmrForLevel(original_file, 1) 

        test_file = "../_output/fort.q" + str(i).zfill(4)
        test_case = ramr.ReadAmrForLevel(test_file, 1)

        error = finderror.error_between_geoclaw(test_case, original_case, "relative")
        error_pct = finderror.error_between_geoclaw(test_case, original_case, "percent")
        #error_water = np.ma.masked_array(error, mask=df_original["total_height"] == 0)
        print "timestep = " + str(i)
        print "L_2 norm = ", norm(error.water[~error.water.mask]) 
        print "L_inf norm = ", norm(error.water[~error.water.mask],np.inf)

        #Plot parameters
        if ploterror:
            plot_title = "Error at timestep " + str(i)
            vmin = -0.0001
            vmax = 0.0001
            savefile = "error_state"+str(i) + ".pdf"
            plotmap.class_contour(error, plot_title, vmin, vmax, savefile=savefile)
        
        print "Max error percentage = ",norm(error_pct.water[~error_pct.water.mask],np.inf)
        print ""
        #plt.figure(i+6)
        if ploterror:
            plot_title = "timestep = " + str(i)
            vmin = -20.0
            vmax = 20.0
            savefile="errorpct_state"+str(i)
            plotmap.class_contour(error_pct, plot_title, vmin, vmax,savefile=None)


if __name__ == "__main__":
    #xlower = -98.0
    #xupper = 98.0
    #ylower = -98.0
    #yupper = 98.0
    #mx = 50
    #my = 50
    ictype = "hump"
    num_time_steps = 10

    #X = np.linspace(xlower, xupper, mx)
    #Y = np.linspace(ylower, yupper, my)
    #xv,yv = np.meshgrid(X,Y)


    error_calc(ictype, num_time_steps, ploterror=True)
    #calcerror(ictype, num_time_steps, ploterror=True)
