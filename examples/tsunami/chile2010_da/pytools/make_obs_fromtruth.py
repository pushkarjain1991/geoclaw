import numpy as np
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import read_amr
import glob

TRUTH_PATH = '../truth/_output/'

def extend_domain(pandas_df, num_cells_extended=200):
    """Extend the domain by num_cells_extended in southwest corner
    """
    xs = pandas_df.drop_duplicates('xcoord').xcoord
    ys = pandas_df.drop_duplicates('ycoord').ycoord
    spacing = xs[1]-xs[0]
    x2 = min(pandas_df.xcoord)
    y2 = min(pandas_df.ycoord)
    x1 = x2 - num_cells_extended*spacing
    y1 = y2 - num_cells_extended*spacing
    x_extended = np.arange(x1, x2, spacing)
    y_extended = np.arange(y1, y2, spacing)
    extended_mesh_grid_x, extended_mesh_grid_y = np.meshgrid(x_extended, y_extended)
    extended_dict = {'xcoord':np.ravel(extended_mesh_grid_x), 'ycoord':np.ravel(extended_mesh_grid_y)}
    extended_df = pd.DataFrame(data=extended_dict)
    pandas_df = pd.concat([extended_df, pandas_df],ignore_index=True)
    return pandas_df


def obs_plot(pandas_df, obs_x, obs_y, vmin=-0.2, vmax=0.2):
    import matplotlib.tri as tr
    from clawpack.visclaw import geoplot
    pcolor_cmap_water = geoplot.tsunami_colormap
    data = np.zeros(len(pandas_df))
    #Mask land
    mask_land = pandas_df.height < 1.0e-6
    mask_land = mask_land.values

    triang = tr.Triangulation(pandas_df.xcoord, pandas_df.ycoord)
    mask = np.all(np.where(mask_land[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)

    fig, ax = plt.subplots(1,1)
    ax.patch.set_color('green')
    plt.tricontourf(triang, data, 40, vmin=vmin, vmax=vmax, cmap=plt.get_cmap(pcolor_cmap_water))
    plt.scatter(obs_x, obs_y, s=5, c='black')
    ax.set_xlim([-130.0, -60.0])
    ax.set_ylim([-70.0, 0.0])
    plt.title('Observation grid')
    plt.savefig('observation.pdf')

def make_obs_from_truth(obs_times, obs_x, obs_y):
    """
    Construct observation from truth.
    User supplies the times and positions where the observations are recorded.
    The observation data is compared to the truth output and obs_step* is created.
    """
    
    #Find all fort.t files
    fortt_file_list = sorted(glob.glob(TRUTH_PATH+"fort.t*"))
    #Find all fort.q files
    fortq_file_list = sorted(glob.glob(TRUTH_PATH+"fort.q*"))

    # Read all the fort.t files and store the time in times_in_file
    times_in_file = []
    for fortt_file in sorted(fortt_file_list):
        times_in_file.append(get_time_in_file(fortt_file))
        times_in_file = sorted(times_in_file)
    
    # Obtain files where user time matches fort.t time
    file_index = np.searchsorted(times_in_file, obs_times)

    # Store all the interested fort.q files in new_file_list
    new_file_list = [fortq_file_list[i] for i in file_index]
    
    for i, fortq_file in enumerate(new_file_list):
        original_case = read_amr.ReadAmr(fortq_file, sort_in_frame=True)
        original_case = original_case.pandas_dataframe
        #original_case = extend_domain(original_case)

        mask1 = np.isclose(original_case['xcoord'].values[:,None], obs_x).any(axis=1)
        mask2 = np.isclose(original_case['ycoord'].values[:,None], obs_y).any(axis=1)
        mask3 = original_case.height > 1.0e-6
        orpd_df = original_case[mask1 & mask2 & mask3]
        
        #Construct observations
        obs_file = "obs_step"+str(i+1)+".txt"
        orpd_df.to_csv(obs_file, sep="\t", header=False, index=False, columns=["xcoord", "ycoord", "eta"])
        print "user time:", obs_times[i], "; fortt time:", times_in_file[file_index[i]], "; Read geoclaw output ", fortq_file, " to produce ", obs_file

    obs_plot(original_case, orpd_df.xcoord.values, orpd_df.ycoord.values)

def get_time_in_file(fortt_file):
    """
    Read the fort.t file and obtain the time
    """
    with open(fortt_file,'r') as f:
        content = f.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    content = [x.strip() for x in content]
    # Obtaining the first string and converting to float
    time_in_file = float(content[0].split()[0])
    return time_in_file


if __name__=="__main__":
    #xobs_start = -119.7
    #xobs_end = -71.7
    #yobs_start = -59.7
    #yobs_end = -19.7
    #nxobs = 21
    #nyobs = 21
    
    #xobs_start = -113.7
    #xobs_end = -83.7
    #yobs_start = -47.7
    #yobs_end = -7.7
    #nxobs = 11
    #nyobs = 11
    
    obs_times = np.arange(10800.,46801,900)

    #obs_x = np.linspace(xobs_start,xobs_end, nxobs)
    #obs_y = np.linspace(yobs_start,yobs_end, nyobs)
    # 52 obs
    #obs_x = np.arange(-113.7, -69.6, 6.0)
    #obs_y = np.arange(-47.7,-7.7, 6.0)
    # 42 obs
    obs_x = np.arange(-113.7, -83.6, 6.0)
    obs_y = np.arange(-47.7,-7.7, 6.0)
    # 35 obs
    #obs_x = np.arange(-113.7, -89.6, 6.0)
    #obs_y = np.arange(-47.7,-7.7, 6.0)
    # 12 obs
    #obs_x = np.arange(-113.7, -83.6, 12.0)
    #obs_y = np.arange(-47.7,-7.7, 12.0)
    #4 obs
    #obs_x = np.arange(-113.7, -89.6, 24.0)
    #obs_y = np.arange(-47.7,-7.7, 30.0)
    
    make_obs_from_truth(obs_times, obs_x, obs_y)
