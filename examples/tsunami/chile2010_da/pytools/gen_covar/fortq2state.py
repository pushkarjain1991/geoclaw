import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
import sys
sys.path.insert(0,'../')
import read_amr2

LARGE_ENS_PATH = "/h2/pkjain/Desktop/Pushkar/clawpack/geoclaw/examples/tsunami/chile2010_da/Pre_ens_gen/_output_100/"
def fortq2state_temporal(ens_num, time_step_interested):
    steps = 400
    for timestep in range(1,steps+1):
        infile = "../../_output/fort.q" + str(timestep).zfill(4)
        #outfile = "../_output/state.q" + str(timestep).zfill(4)
        outfile = "../../_output/state.q" + str(timestep)
        read_class = read_amr2.ReadAmr2(infile, sort_in_frame=False)
        dataframe = read_class.pandas_dataframe
        print "timestep = ", timestep
        dataframe.eta.to_csv(outfile, index=False)

def fortq2state_state(ens_num, time_step_interested):
    infile = LARGE_ENS_PATH + "_output_" + str(ens_num) + "_for/fort.q" + str(time_step_interested).zfill(4)
    outfile = LARGE_ENS_PATH + "state.q" + str(ens_num)
    read_class = read_amr2.ReadAmr2(infile, sort_in_frame=False)
    dataframe = read_class.pandas_dataframe
    print "Ens = ", ens_num

    #Mask land region. If not done, perts will be added in fault
    dataframe.ix[np.isclose(dataframe.height.values,0.0), 'eta'] = 0.0
    dataframe.eta.to_csv(outfile, index=False)
    return dataframe.eta.values

def plot_covar(H):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('colorMap')
    plt.imshow(H)
    ax.set_aspect('equal')

    cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
    cax.get_xaxis().set_visible(False)
    cax.get_yaxis().set_visible(False)
    cax.patch.set_alpha(0)
    cax.set_frame_on(False)
    plt.colorbar(orientation='vertical')
    plt.savefig('covar_mat')


if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 1):
        sys.stdout.write('No command line dim_ens given. Using 100')
        dim_ens = 100
    else:
        dim_ens = int(sys.argv[1])
        sys.stdout.write('Using dim_ens ' + str(dim_ens))
    
    time_step_interested = 8

#    for ens_num in range(dim_ens):
#        fortq2state('space', ens_num, time_step_interested)
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores)(delayed(fortq2state_state)(ens_num, time_step_interested) for ens_num in range(dim_ens))
