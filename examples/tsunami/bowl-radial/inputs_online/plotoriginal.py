import plotmap
import ReadAmrForLevel as ramr

if __name__=="__main__":
    obs_time_list = [20,40,60,80,100,120,140,160,180]

    for i,j in enumerate(obs_time_list):
        read_geoclaw_output = "../_output_original_hump/fort.q" + str(i+1).zfill(4)
        original_case = ramr.ReadAmrForLevel(read_geoclaw_output, 1.0)
        plot_title = "Actual state: Timestep " + str(j)
        vmin = -0.5
        vmax = 0.5
        savefile = "original_state"+str(i+1)+".pdf"
        print "Saving", savefile
        plotmap.class_contour(original_case,plot_title, vmin,vmax,savefile=savefile)
        

