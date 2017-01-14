import numpy as np
from joblib import Parallel, delayed
import multiprocessing
import grid_plot as grplot


def init_ens_plot(ens_num, init_ens, plot_var, fortfile = '../_output/fort.q0005'):
    import matplotlib.pyplot as plt
    import matplotlib.tri as tr
    from read_amr2 import ReadAmr2
    from clawpack.visclaw import geoplot
    pcolor_cmap_water = geoplot.tsunami_colormap

    file1 = ReadAmr2(fortfile, sort_in_frame=False)
    file1 = file1.pandas_dataframe


    #Mask land
    #file1.ix[file1.height==0.0, 'eta'] = np.nan
    mask_land = file1.height < 1.0e-6
    mask_land = mask_land.values
    triang = tr.Triangulation(file1.xcoord, file1.ycoord)
    mask = np.all(np.where(mask_land[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)

    fig, ax = plt.subplots(1,1)
    plt.gca().patch.set_color('green')
    plt.tricontourf(triang, init_ens, 20, cmap = plt.get_cmap(pcolor_cmap_water),vmin = -0.2, vmax = 0.2)
    plt.colorbar()
    plt.title('Initial' + plot_var)
    plt.savefig('init_' + plot_var + str(ens_num) + '.png')
    print "Ens num", ens_num


if __name__ == "__main__":
    # *** Inputs ***
    dim_ens = 9
    init_pert_array = np.loadtxt("../_output/check_obs.txt")
    init_ens_array = [np.loadtxt("../_output/init_ens" + str(i)) for i in range(1,dim_ens+1)]

    ens_num_array = np.arange(dim_ens)

    # *** Plot single member perts***
    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(init_ens_plot)(i,j,'pert') for i,j in zip(ens_num_array, init_pert_array))


    # *** Construct pert grid plot ***
    image_array = ['init_pert' + str(ens_num) + '.png' for ens_num in np.arange(dim_ens)]
    outfile = 'init_pert_grid.png'
    grplot.grid_plot(image_array, outfile)
    print 'Saved init_pert_grid.png'
    
#    # *** Plot single member state***
#    num_cores = multiprocessing.cpu_count()
#    results = Parallel(n_jobs=num_cores)(delayed(init_ens_plot)(i,j,'state') for i,j in zip(ens_num_array, init_ens_array))
#    
#    # *** Construct state grid plot ***
#    image_array = ['init_state' + str(ens_num) + '.png' for ens_num in np.arange(dim_ens)]
#    outfile = 'init_state_grid.png'
#    grplot.grid_plot(image_array, outfile)
#    print 'Saved init_state_grid.png'
