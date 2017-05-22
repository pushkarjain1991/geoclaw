import matplotlib
matplotlib.use('Agg')
#import matplotlib as mpl
#mpl.use('pgf')
#pgf_with_latex = {                      # setup matplotlib to use latex for output
#    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
#    "text.usetex": True,                # use LaTeX to write all text
#    "font.family": "serif",
#    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
#    "font.sans-serif": [],
#    "font.monospace": [],
#    "axes.labelsize": 10,               # LaTeX default is 10pt font.
#    "text.fontsize": 10,
#    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
#    "xtick.labelsize": 8,
#    "ytick.labelsize": 8,
#    "pgf.preamble": [
#        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
#        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
#        ]
#    }
#mpl.rcParams.update(pgf_with_latex)
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#from matplotlib2tikz import save as tikz_save

singular_value_file = 'singular_values'
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
plt.savefig('Eigen.pdf')
#plt.savefig('Eigen.pgf')

