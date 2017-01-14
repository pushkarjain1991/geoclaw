import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

singular_value_file = '../singular_values'
num_eval_plot = 12

#*** Read the eigen values ***
eigen_values = np.loadtxt(singular_value_file)

#*** Caluclate explained variance ***
tot = sum(eigen_values)
var_exp = [(i / tot)*100 for i in eigen_values]
cum_var_exp = np.cumsum(var_exp)


#*** Plotting ***
def autolabel(rects, ax):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height), ha='center', va='bottom')


fig, ax = plt.subplots(nrows=2,ncols=1)
plt.suptitle('Eigen value distribution', fontsize=14, fontweight='bold')

# Plot 1 - Eigen values
ax[0].plot(eigen_values)
ax[0].set_ylim(bottom=0.0)
ax[0].set_ylabel('Eigen value')

#Plot explained variance
rects = ax[1].bar(np.arange(num_eval_plot), var_exp[:num_eval_plot])
autolabel(rects, ax=ax[1])
width = rects[0].get_width()/2.0
ax[1].plot(np.arange(num_eval_plot) + 0.5, cum_var_exp[:num_eval_plot],'-ro')

ax[1].set_xticks(np.arange(num_eval_plot)+width)
ax[1].set_xticklabels(np.arange(num_eval_plot))
ax[1].set_ylabel('% Explained variance')

ax[1].set_xlabel('Eigen value index')
plt.show()

