import ReadAmrForLevel as ramr
import read_amr2
#import finderror
import matplotlib.pyplot as plt
import numpy as np

def error_calc(ictype, num_time_steps):
    print "Domain is reduced"
    print "Timestep, L2 norm, L_inf norm"
    for i in range(num_time_steps):
        #original_file = "../_output_original_" + ictype + "/fort.q" + str(i).zfill(4)
        original_file = "../_output_original/fort.q" + str(i).zfill(4)
        original_case = read_amr2.ReadAmr2(original_file)
        original_case = original_case.amrdataframe()
        original_case = original_case[np.isclose(original_case.amrlevel,1)]

        test_file = "../_output/fort.q" + str(i).zfill(4)
        test_case = read_amr2.ReadAmr2(test_file)
        test_case = test_case.amrdataframe()
        test_case = test_case[np.isclose(test_case.amrlevel,1)]
        #test_case = test_case[test_case.amrlevel == 1]
        #test_case = ramr.ReadAmrForLevel(test_file, 1)


        #mask1 = original_case["xcoord"].isin(test_case.xcoord)
        #mask2 = original_case["ycoord"].isin(test_case.ycoord)
        mask1 = np.isclose(original_case['xcoord'].values[:,None], np.unique(test_case['xcoord'].values)).any(axis=1)
        mask2 = np.isclose(original_case['ycoord'].values[:,None], np.unique(test_case['ycoord'].values)).any(axis=1)
        subset_case = original_case[mask1 & mask2]
        subset_case = subset_case.reset_index(drop=True)


        subset_case = subset_case[(subset_case.xcoord <= -85.0)]
        test_case = test_case[(test_case.xcoord <= -85.0)]

        assert len(subset_case) == len(test_case)

        # Error analysis
        error_eta = subset_case.height - test_case.height
        plt.tricontourf(subset_case.xcoord, subset_case.ycoord, error_eta, 20)
        plt.colorbar()
        plt.show()
        print str(i), np.linalg.norm(error_eta), np.linalg.norm(error_eta, np.inf)

if __name__ == "__main__":
    ictype = ""
    num_time_steps = 13

    #X = np.linspace(xlower, xupper, mx)
    #Y = np.linspace(ylower, yupper, my)
    #xv,yv = np.meshgrid(X,Y)


    error_calc(ictype, num_time_steps)
    #calcerror(ictype, num_time_steps, ploterror=True)
