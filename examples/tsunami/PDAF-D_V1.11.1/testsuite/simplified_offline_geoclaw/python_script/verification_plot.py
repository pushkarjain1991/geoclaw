import ReadAmrForLevel as ramr
# import finderror as errg
import plotmap
import numpy as np
import finderror


def verify_das(dtobs, ens_number, output_times, num_ens):
    for k, j in enumerate(dtobs[:-1]):

        # Set the test_case file to be compared to original geoclaw
        read_verification_output = "./ens_" + str(ens_number) + "_" + str(k) \
                                    + "/fort.q00" + str(output_times)
        test_case = ramr.ReadAmrForLevel(read_verification_output, 1.0)

        # Set the original geoclaw
        original_fortq_time = str((np.size(dtobs)-1) * output_times *
                                  int(dtobs[k+1])/int(dtobs[-1])).zfill(2)
        original_fortq_file = "./original_geoclaw/fort.q00"+original_fortq_time
        print "\nOriginal file read is " + original_fortq_file
        original_case = ramr.ReadAmrForLevel(original_fortq_file, 1.0)

        error_class = finderror.error_between_geoclaw(test_case,
                                                      original_case,
                                                      type="relative")
        error_percent_class = finderror.error_between_geoclaw(test_case,
                                                              original_case,
                                                              type="percent")

        plotmap.class_contour(original_case, "Original WSE: ens_number = " +
                              str(ens_number) + "; num_ens = " + str(num_ens) +
                              "; time = " + str(dtobs[k+1]), -0.9, 0.9)

        plotmap.class_contour(test_case, "WSE: ens_number = " +
                              str(ens_number) + "; num_ens = " +
                              str(num_ens) + "; time = " + str(dtobs[k+1]),
                              -0.9, 0.9)
        plotmap.class_contour(error_class, "WSE error: ens_number = " +
                              str(ens_number) + "; num_ens = " + str(num_ens) +
                              "; time = " + str(dtobs[k+1]), -0.01, 0.01)
        plotmap.class_contour(error_percent_class, "WSE error %: \
                              ens_number = " + str(ens_number)+"; num_ens = \
                              " + str(num_ens) + "; time = " +
                              str(dtobs[k+1]), -100.0, 100.0)
 #       plotmap.docontour(test_case.mxv, test_case.myv, error_reshaped_masked_eta_water, error_reshaped_masked_eta_land, "WSE error: ens_number = "+str(ens_number)+"; num_ens = " + str(num_ens)+"; time = " + str(dtobs[k+1]), -0.9,0.9)
  #      plotmap.docontour(test_case.mxv,test_case.myv,error_reshaped_masked_eta_water_percent, error_reshaped_masked_eta_land, "WSE error %: ens_number = "+str(ens_number)+"; num_ens = " + str(num_ens)+"; time = " + str(dtobs[k+1]), -200.0,200.0)
