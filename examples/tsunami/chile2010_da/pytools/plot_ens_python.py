from PIL import Image
import subprocess
from joblib import Parallel, delayed
import multiprocessing
import grid_plot as grplot
import numpy as np

def processInput(folder_num):
    """
    Construct plot for a single ensemble member by calling the shell subprocess
    """
    filename = '../_output/_output_' + str(folder_num) + '_for/'
    shell_str = 'python /workspace/Pushkar/clawpack/visclaw/src/python/visclaw/plotclaw.py ../_output/_output_' + str(folder_num) + '_for _plots_' + str(folder_num) + ' ../setplot.py'
    subprocess.call(shell_str, shell=True)

if __name__ == "__main__":
    dim_ens = 9
    timesteps = np.arange(12,53,1)
    
    # Construct plot for every ensemble
    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in range(dim_ens))

    #Construct grid for every time step
    for t in timesteps:
        image_array = ['_plots_' + str(ens) + '/frame' + str(t).zfill(4) +'fig0.png' for ens in range(dim_ens)]
        outfile = 'grid_img_' + str(t) + '.png'
        grplot.grid_plot(image_array, outfile)

    # Construct animated gif
    print "Making Gif..."
    image_array =  ['grid_img_' + str(t) + '.png' for t in timesteps]
    outfile = 'ens_grid.gif'
    grplot.grid2gif2(image_array, outfile)

    print "----DONE----"
