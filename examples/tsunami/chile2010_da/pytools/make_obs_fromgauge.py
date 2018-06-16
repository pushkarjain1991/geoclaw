import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

x_coord = -86.392
y_coord = -17.975
#bath = 4374.368
bath = 0.0
dt = 60.0
tstart = 10800.0
tstop = 46800.0
#nstops = 21
#tstart = 1800.0
#tstop = 46800.0
#nstops = 26
obs_data = pd.read_csv('../32412_notide.txt', names=['time','wse'], index_col=False, sep=" " )

#obs_time_list = np.linspace(tstart, tstop, nstops)
obs_time_list = np.arange(10800.0, 46801.0, 900.0)

extracted_data = obs_data[np.isclose(obs_data['time'].values[:,None], obs_time_list).any(axis=1)]

first_val = extracted_data.groupby('time', as_index=False).first()
first_val['timestep'] = first_val.time/3600.0
print first_val
#first_val.plot('time','wse',ax=ax1)

mean_val = extracted_data.groupby('time', as_index=False).mean()

obs_timestep_list = first_val.time/dt

#fig, ax1 = plt.subplots()
#mean_val.plot('time', 'wse',ax=ax1)
#plt.legend(['first', 'mean'])
#plt.show()

wse_array = first_val['wse'].values + bath

for i,j in enumerate(wse_array):
    #savefile = "obs_step" + str(int(obs_timestep_list[i])) + ".txt"
    savefile = "obs_step" + str(i+1) + ".txt"
    with open(savefile, 'w') as f1:
        f1.write(str(x_coord) + " " + str(y_coord) + " " + str(j))
        #f1.write(str(x_coord) + " " + str(y_coord) + " " + str(2.0))
        f1.write("\n")
        #f1.write("-110.0 -40.0 .01")

#for i,j in enumerate(obs_time_list):
#    obs_file = "obs_step"+str(int(j))+".txt"
    
    
