import numpy as np
import pandas as pd
from read_amr import ReadAmr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tr
from joblib import Parallel, delayed
import multiprocessing
import grid_plot as grplot
from PIL import Image
import general_plot as gp

def var_plot(varfile, timestep, fortfile = '../_output/fort.q0005'):
    
    print "Plotting stddev at timestep ", timestep
    fortfile = "../Assimilated_results/localization_effect/35_obs/radius_15/_output/fort.q0005"

    file1 = ReadAmr(fortfile, sort_in_frame=False)
    file1 = file1.pandas_dataframe

    var = np.loadtxt(varfile)
    stddev = np.sqrt(var)

    title = 'Standard deviation | Assimilation step ' + str(timestep)
    savefile = 'var' + str(timestep).zfill(4) + '.pdf'
    gp.general_plot(file1, stddev, 0.0, 0.03, title = title, savefile = savefile, colorbar=True)


def plot_parallel(timesteps, var_files):
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores)(delayed(var_plot)(varfile,timestep) for varfile,timestep in zip(var_files, timesteps))


if __name__ == "__main__":
    import imageio
    #timesteps = np.arange(3, 12)
    timesteps = np.arange(1, 41)
    
    # *** Plot stddev ***
    #var_files = ["../_output/variance_" + str(i) + "_ana.txt"  for i in timesteps]
    var_files = ["../Assimilated_results/localization_effect/35_obs/radius_15/_output/variance_" + str(i) + "_ana.txt"  for i in timesteps]
    plot_parallel(timesteps, var_files)

    # *** Construct GIF ***
    print "Making GIF... "
    #var_images_array = []
    #var_images_list = ["var" + str(i) + ".png" for i in timesteps]
    #for var_image in var_images_list:
    #    im = imageio.imread(var_image)
    #    var_images_array.append(im)
    grplot.grid2gif4("var*.pdf", 'variance.gif')

    print "Process done"
