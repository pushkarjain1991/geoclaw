import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tr
from clawpack.visclaw import geoplot
import matplotlib.colors as mcol


def general_plot(pandas_df, data, vmin, vmax, title=None, savefile=None, ax=None, colorbar=False):
    pcolor_cmap_water = geoplot.tsunami_colormap

    #Mask land
    #file1.ix[file1.height==0.0, 'eta'] = np.nan
    mask_land = pandas_df.height < 1.0e-6
    mask_land = mask_land.values
    triang = tr.Triangulation(pandas_df.xcoord, pandas_df.ycoord)
    mask = np.all(np.where(mask_land[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)

    if ax is None:
        fig, ax = plt.subplots(1,1)
    ax.patch.set_color('green')
    ax = ax.tricontourf(triang, data, 100, vmin=vmin, vmax=vmax, cmap = plt.get_cmap(pcolor_cmap_water))
    if colorbar:
        cbar = fig.colorbar(ax)
    if title is not None:
        plt.title(title)
    if savefile is not None:
        plt.savefig(savefile)

    return ax
