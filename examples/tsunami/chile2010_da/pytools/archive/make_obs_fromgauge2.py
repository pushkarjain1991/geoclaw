import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

x_coord = -86.392
y_coord = -17.975
bath = 4374.368
obs_data = pd.read_csv('32412_notide.txt', names=['timestep','wse'], index_col=False, sep=" " )

obs_time_list = np.linspace(30.0, 780.0, 26)

print obs_data
extracted_data = obs_data[np.isclose(obs_data['time'].values[:,None], obs_time_list).any(axis=1)]
print extracted_data

fig, ax1 = plt.subplots()
first_val = extracted_data.groupby('time').first()
first_val.plot(ax=ax1)

mean_val = extracted_data.groupby('time').mean()
mean_val.plot(ax=ax1)

plt.legend(['first', 'mean'])
plt.show()

wse_array = first_val['wse'].values + bath

for i,j in enumerate(wse_array):
    savefile = "obs_step" + str(int(obs_time_list[i])) + ".txt"
    with open(savefile, 'w') as f1:
        f1.write(str(x_coord) + " " + str(y_coord) + " " + str(j))
        f1.write("\n")
        #f1.write("-110.0 -40.0 3840.0")

#for i,j in enumerate(obs_time_list):
#    obs_file = "obs_step"+str(int(j))+".txt"
    
    
