import pandas as pd
from read_amr import ReadAmr
from general_plot import general_plot
import numpy as np
import matplotlib.pyplot as plt
from grid_plot import grid2gif4

def plot_err(dir1, timestep1, dir2, timestep2, plotting=True):
    """
    Error between mean state and truth
    """
    #First directory
    file1 = dir1 + "fort.q" + str(timestep1).zfill(4)
    df1 = ReadAmr(file1, sort_in_frame=True)
    df1 = df1.pandas_dataframe
    df1 = df1[np.isclose(df1.amrlevel,1.0)]
    
    #Second directory
    file2 = dir2 + "fort.q" + str(timestep2).zfill(4)
    df2 = ReadAmr(file2, sort_in_frame=True)
    df2 = df2.pandas_dataframe
    df2 = df2[np.isclose(df2.amrlevel,1.0)]
    df2.sort_values(by=['ycoord','xcoord'], ascending=[True, True], inplace=True)
    
    if 0:
        mask1 = (df1.xcoord > -76) & (df1.xcoord < -68) & (df1.ycoord > -40) & (df1.ycoord < -32)
        mask2 = (df2.xcoord > -76) & (df2.xcoord < -68) & (df2.ycoord > -40) & (df2.ycoord < -32)
        df1 = df1[~mask1]
        df2 = df2[~mask2]
    
    #Error
    error = df1.eta.values - df2.eta.values
    error_norm = np.linalg.norm(error)
    
    #Plot 
    if plotting:
        plot_title = "Error plot at assimilation step = " + str(timestep1)
        savefile = "plot_err" + str(timestep1) + ".pdf"
        general_plot(df1, error, vmin=-0.25, vmax=0.25, title=plot_title, savefile=savefile)
        print "Generated " + savefile
    
    return error_norm


if __name__=="__main__":

    # Truth
    dir1 = "../truth/_output/"
    timesteps1 = np.arange(12,53,1)

    # Assimilated state
    dir2 = "../_output/"
    #dir2 = "/h2/pkjain/Desktop/Pushkar/clawpack/geoclaw/examples/tsunami/chile2010_da/Assimilated_results/localization_effect/4_obs/radius_"+ str(i) + "/_output/"
    timesteps2 = np.arange(8,49,1)
    
    # Free run state
    dir3 = "../free_run/_output/"
    timesteps3 = np.arange(12,53,1)

    # Calculate error of assimilated state
    error_array2 = []
    for counter, (timestep1, timestep2) in enumerate(zip(timesteps1, timesteps2)):
        err = plot_err(dir1, timestep1, dir2, timestep2, plotting=True)
        error_array2.append(err)

    # Calculate error of free run
    error_array1 = []
    for counter, (timestep1, timestep2) in enumerate(zip(timesteps1, timesteps3)):
        err = plot_err(dir1, timestep1, dir3, timestep2, plotting=False)
        error_array1.append(err)

    fig,ax = plt.subplots(1,1)
    p1 = ax.plot(timesteps1, error_array2, '-bo', label='Assimilation run')
    p2 = ax.plot(timesteps1, error_array1, '-ro', label='Free run')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Error norm')
    plt.title('Error comparison for twin experiment')
    plt.legend(loc='best')
    plt.savefig('error_check_enkf.pdf')
    print "Plotted error norm plot error_check_enkf.pdf"

    #filename = 'reduced_error_local' + str(i) + '.txt'
    filename = 'reduced_error_local.txt'
    with open(filename, 'w') as f1:
        for item in error_array2:
            print>>f1, item

    with open('error_freerun.txt', 'w') as f1:
        for item in error_array1:
            print>>f1, item
    
    print "Making GIF... "
    grid2gif4("plot_err*.pdf", 'plot_err.gif')

