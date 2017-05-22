import matplotlib as mpl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.tri as tr
import matplotlib.colors as mcol
import general_plot as gp
from read_amr import ReadAmr
from clawpack.visclaw import geoplot

# input
dim_ens = 9
init_pert_array = np.loadtxt("../_output/check_obs.txt")

fortfile = '../_output/fort.q0000'
fortdata = ReadAmr(fortfile, sort_in_frame=False)
pandas_df = fortdata.pandas_dataframe

vmin = min(np.amin(init_pert_array, axis=1))
vmax = max(np.amax(init_pert_array, axis=1))
norm_val = max(abs(vmin), abs(vmax))
if (vmin < 0) and (vmax > 0):
    vmin = -1*norm_val
    vmax = norm_val
else: 
    print "Nothing happened"
vmin = -0.3
vmax = 0.3

fig, axes = plt.subplots(nrows=3,ncols=3, figsize=(15,15))
flat_axes = axes.flatten
for ax, ens_num, data in zip(axes.flat, np.arange(dim_ens), init_pert_array):
    ax = gp.general_plot(pandas_df, data, vmin, vmax, ax=ax)

# Adding common colorbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
cmap = plt.get_cmap(geoplot.tsunami_colormap)
norm = mcol.Normalize(vmin=vmin,vmax=vmax)
#bounds = np.linspace(vmin,vmax,50)
ticks = np.linspace(vmin,vmax,11)
cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, ticks=ticks, norm=norm)
cb1.set_label('Water surface elevation perturbation', size=26)
cb1.ax.tick_params(labelsize=22) 

plt.savefig('Initial_ensemble.pdf')
print "Plotted initial ensemble in Initial_ensemble.pdf"

