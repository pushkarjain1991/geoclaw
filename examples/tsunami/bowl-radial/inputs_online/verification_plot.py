import ReadAmrForLevel as ramr
# import finderror as errg
import matplotlib.pyplot as plt
import plotmap
import numpy as np
import finderror


def verify_das(dtobs, num_ens):
    for k, j in enumerate(range(dtobs)):

        # Set the test_case file to be compared to original geoclaw
        read_verification_output = "../_output/fort.q000" + str(j)
        test_case = ramr.ReadAmrForLevel(read_verification_output, 1.0)

        # Set the original geoclaw
        original_fortq_file = "../_output_original_hump/fort.q000" + str(j)

        print "Original file read is " + original_fortq_file
        original_case = ramr.ReadAmrForLevel(original_fortq_file, 1.0)

        error_class = finderror.error_between_geoclaw(test_case,
                                                      original_case,
                                                      type="relative")
        #np.savetxt("hum1", error_class.get_water())
        error_percent_class = finderror.error_between_geoclaw(test_case,
                                                              original_case,
                                                              type="percent")

        plotmap.class_contour(original_case, "Original WSE ; num_ens = " + str(num_ens) +
                              "; time = " + str(j), -0.9, 0.9)

        plotmap.class_contour(test_case, "WSE ; num_ens = " +
                              str(num_ens) + "; time = " + str(j),
                              -0.9, 0.9)
        plotmap.class_contour(error_class, "WSE error ; num_ens = " + str(num_ens) +
                              "; time = " + str(j), -0.01, 0.01)
        plotmap.class_contour(error_percent_class, "WSE error %: num_ens = " + str(num_ens) + "; time = " +  str(j), -100.0, 100.0)
        plt.show()
 #       plotmap.docontour(test_case.mxv, test_case.myv, error_reshaped_masked_eta_water, error_reshaped_masked_eta_land, "WSE error: ens_number = "+str(ens_number)+"; num_ens = " + str(num_ens)+"; time = " + str(dtobs[k+1]), -0.9,0.9)
  #      plotmap.docontour(test_case.mxv,test_case.myv,error_reshaped_masked_eta_water_percent, error_reshaped_masked_eta_land, "WSE error %: ens_number = "+str(ens_number)+"; num_ens = " + str(num_ens)+"; time = " + str(dtobs[k+1]), -200.0,200.0)

if __name__ == "__main__":
    verify_das(5, 5)
