import numpy as np
import pandas as pd
from read_amr2 import ReadAmr2
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tr
from joblib import Parallel, delayed
import multiprocessing
import grid_plot as grplot

def var_plot(varfile, timestep, fortfile = '../_output/fort.q0005'):
    
    print "Plotting stddev at timestep ", timestep

    file1 = ReadAmr2(fortfile, sort_in_frame=False)
    file1 = file1.pandas_dataframe

    var = np.loadtxt(varfile)
    stddev = np.sqrt(var)

    #Mask land
    #file1.ix[file1.height==0.0, 'eta'] = np.nan
    mask_land = file1.height < 1.0e-6
    mask_land = mask_land.values
    triang = tr.Triangulation(file1.xcoord, file1.ycoord)
    mask = np.all(np.where(mask_land[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)

    fig, ax = plt.subplots(1,1)
    plt.gca().patch.set_color('green')
    plt.tricontourf(triang, stddev, 20, vmin = 0.0, vmax = 0.01)
    plt.colorbar()
    plt.title('Standarad deviation | Timestep ' + str(timestep))
    plt.savefig('var' + str(timestep) + '.png')

def plot_serial():
    timesteps = np.arange(3, 22)
    var_files = ["../_output/variance_" + str(i) for i in timesteps]
    timesteps = np.arange(3, 23)
    for var_file, timestep in zip(var_files, timesteps):
        var_plot(var_file, timestep)

def plot_parallel(timesteps, var_files):
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores)(delayed(var_plot)(i,j) for i,j in zip(var_files, timesteps))


if __name__ == "__main__":
    timesteps = np.arange(3, 12)
    
    # *** Plot stddev ***
    var_files = ["../_output/variance_" + str(i) for i in timesteps]
    plot_parallel(timesteps, var_files)

    # *** Construct GIF ***
    print "Making GIF... "
    var_images_array = ["var" + str(i) + ".png" for i in timesteps]
    grplot.grid2gif2(var_images_array, 'variance.gif')

    print "Process done"
