import ensemble_class
import subprocess
import os

def run_geoclaw(k, dtobs, mx, my, geoclaw_input, topography, output_times, max_amr, geoclaw_exec):

    hello = ensemble_class.ensemble()
    hello.rundata.clawdata.t0 = dtobs[k]
    hello.rundata.clawdata.tfinal = dtobs[k+1]
    hello.rundata.clawdata.num_cells[0] = mx
    hello.rundata.clawdata.num_cells[1] = my
    print mx, my
    hello.rundata.qinit_data.qinitfiles[-1][-1]=geoclaw_input
    #hello.rundata.qinit_data.qinitfiles[-1]=[1,2,geoclaw_input]
    #hello.rundata.topo_data.topofiles[-1]=[2, 1, 1, 0., 1.e10, topo_path]
    hello.rundata.topo_data.topofiles[-1][-1]=topography
    hello.rundata.clawdata.num_output_times = output_times
    hello.rundata.amrdata.amr_levels_max=max_amr
    hello.rundata.write()
    FNULL = open(os.devnull,'w')
    subprocess.call(geoclaw_exec)
    #subprocess.call(geoclaw_exec, stdout=FNULL, stderr=subprocess.STDOUT)
