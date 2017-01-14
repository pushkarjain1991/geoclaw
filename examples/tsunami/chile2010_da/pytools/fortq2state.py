import read_amr2
import numpy as np
import matplotlib.pyplot as plt


def fortq2state(which_type, ens_num, time_step_interested):
    if which_type is "temporal":
        steps = 400
        for timestep in range(1,steps+1):
            infile = "../_output/fort.q" + str(timestep).zfill(4)
            #outfile = "../_output/state.q" + str(timestep).zfill(4)
            outfile = "../_output/state.q" + str(timestep)
            read_class = read_amr2.ReadAmr2(infile, sort_in_frame=False)
            dataframe = read_class.pandas_dataframe
            print "timestep = ", timestep
            dataframe.eta.to_csv(outfile, index=False)
    elif which_type is "space":
        infile = "../_output_ens_gen/_output_" + str(ens_num) + "_for/fort.q" + str(time_step_interested).zfill(4)
        outfile = "../_output_ens_gen/state.q" + str(ens_num)
        read_class = read_amr2.ReadAmr2(infile, sort_in_frame=False)
        dataframe = read_class.pandas_dataframe
        print "Ens = ", ens_num

        #Mask land region. If not done, perts will be added in fault
        dataframe.ix[np.isclose(dataframe.height.values,0.0), 'eta'] = 0.0
        dataframe.eta.to_csv(outfile, index=False)
        return dataframe.eta.values
    else:
        print "Invalid type"

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
    dim_ens = 64
    time_step_interested = 4
    state_size = 10000

    state_matrix = np.empty([state_size, dim_ens])
    for ens_num in range(dim_ens):
        state_matrix[:,ens_num] = fortq2state('space', ens_num, time_step_interested)

    H = np.cov(state_matrix)
    H[np.abs(H) < 1.0e-6] = 0.0
    print "shape of covariance matrix = ", np.shape(H)
    print H.max()
    #plot_covar(H)


    
