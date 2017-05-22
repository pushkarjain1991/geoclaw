import matplotlib
matplotlib.use('Agg')
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
import grid_plot as grplot
from general_plot import general_plot


def init_ens_plot2(ens_num, init_ens, plot_var, fortfile = '../_output/fort.q0005'):
    from read_amr2 import ReadAmr2

    fortdata = ReadAmr2(fortfile, sort_in_frame=False)
    fortdata = fortdata.pandas_dataframe

    savefile = 'init_' + plot_var + str(ens_num) + '.pdf'
    title = 'Initial' + plot_var
    general_plot(fortdata, init_ens, vmin=-0.2, vmax=0.2, title=title, savefile = savefile)
    print "Ens num", ens_num

if __name__ == "__main__":
    # *** Inputs ***
    dim_ens = 9
    init_pert_array = np.loadtxt("../_output/check_obs.txt")
    init_ens_array = [np.loadtxt("../_output/init_ens" + str(i)) for i in range(1,dim_ens+1)]

    ens_num_array = np.arange(dim_ens)

    # *** Plot single member perts***
    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(init_ens_plot2)(i,j,'pert') for i,j in zip(ens_num_array, init_pert_array))

    # *** Construct pert grid plot ***
    image_array = ['init_pert' + str(ens_num) + '.pdf' for ens_num in np.arange(dim_ens)]
    outfile = 'init_pert_grid.pdf'
    grplot.grid_plot(image_array, outfile)
    print 'Saved init_pert_grid.pdf'
    
    # *** Plot single member state***
    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(init_ens_plot2)(i,j,'state') for i,j in zip(ens_num_array, init_ens_array))
    
    # *** Construct state grid plot ***
    image_array = ['init_state' + str(ens_num) + '.pdf' for ens_num in np.arange(dim_ens)]
    outfile = 'init_state_grid.pdf'
    grplot.grid_plot(image_array, outfile)
    print 'Saved init_state_grid.pdf'
