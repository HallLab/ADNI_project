import clean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

savepath = '../results/plots/'
respath  = '../results/'

#### Figure 1
## PLSDA 
qtpad = clean.QT_pad()
qtpad.PLS_DA()
fig = plt.figure(figsize = (12,6))
fig.suptitle('PLS-DA', fontsize=16)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax1.set_xlabel('Component 1', fontsize = 14)
ax1.set_ylabel('Component 2', fontsize = 14)
ax1.set_title('Scores', fontsize = 15)
targets = ['CN',
           'MCI',
           'AD']
color_pastel = cm.get_cmap('Set2')
c = 0
for target in targets:
    indicesToKeep = qtpad.data['DX.bl'] == target
    ax1.scatter(qtpad.scores.loc[indicesToKeep,
                                'Component 1'],
               qtpad.scores.loc[indicesToKeep,
                                'Component 2'],
               color = color_pastel(c),
               s = 50,
               alpha = 0.7)
    c = c + 1
ax1.legend(targets)

ax2.set_xlabel('Component 1', fontsize = 14)
ax2.set_ylabel('Component 2', fontsize = 14)
ax2.set_title('Weights', fontsize = 15)
ax2.scatter(qtpad.x_weights['Weight 1'],
            qtpad.x_weights['Weight 2'])
for segment in qtpad.x_weights.index:
    if segment == 'Ventricles':
        xy_offset = (-10, 10)
    elif segment == 'Entorhinal':
        xy_offset = (-55, -2)
    else:
        xy_offset = (-55, 4)
    ax2.annotate(segment, 
                 xy=(qtpad.x_weights.loc[segment,
                                        'Weight 1'],
                     qtpad.x_weights.loc[segment,
                                        'Weight 2']),
                 xytext=xy_offset,
                 textcoords="offset points")
ax2.scatter(qtpad.y_weights['Weight 1'],
            qtpad.y_weights['Weight 2'])
for group in qtpad.y_weights.index:
    xy_offset = (-10, 5)
    ax2.annotate(group, 
                 xy=(qtpad.y_weights.loc[group,
                                        'Weight 1'],
                     qtpad.y_weights.loc[group,
                                        'Weight 2']),
                 xytext=xy_offset,
                 textcoords="offset points")

filename = savepath +\
           'Figure1.pdf'
plt.savefig(filename,
            dpi=300)

#### Figure 3
# Betas and p-values in NMR
modules_to_highlight = ['green',
                        'turquoise']
met_to_module_names = ['module_colors_nmr.csv',
                       'module_colors_p180.csv']
results_metabolites_nmr = ['results_nmr_Component 1.csv',
                          'results_nmr_Component 2.csv']
results_metabolites_p180 = ['results_p180_Component 1.csv',
                            'results_p180_Component 2.csv']
axis_names = ['Effect in females',
              'Effect in males']
ax_titles = ['Differences in Component 1',
             'Differences in Component 2']

fig = plt.figure(figsize = (12,6))
fig.suptitle('Effect size differences NMR', fontsize=16)
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
axes = [ax1,
        ax2]
for i, name in enumerate(results_metabolites_nmr):
    filename = respath + \
               name
    res = pd.read_csv(filename).\
             set_index(['Variable',
                        'Outcome'])
    mod_name = respath + \
               met_to_module_names[0]
    mods = pd.read_csv(mod_name, 
                       header=None)
    max_value = max(abs(res.loc[:,['Beta_female',
                                   'Beta_male'] ]).max())
    pvals = -np.log(res['pvalue_diff'])
    max_value_round = round(max_value + 0.1, 2)
    min_max   = (-max_value_round,
                 max_value_round)
    pvals = -np.log(res['pvalue_diff'])
    axes[i].scatter(res['Beta_female'],
                    res['Beta_male'],
                    color = 'grey',
                    alpha = 0.4,
                    s = pvals*20)
    for color in modules_to_highlight:
        module_bool = mods[1] == color
        module_metabolites = mods[0][module_bool]
        axes[i].scatter(res.loc[module_metabolites,'Beta_female'],
                        res.loc[module_metabolites,'Beta_male'],
                        color = color,
                        s = pvals[module_metabolites]*20,
                        alpha = 0.8)
    axes[i].set_xlabel(axis_names[0], fontsize = 14)
    axes[i].set_ylabel(axis_names[1], fontsize = 14)
    axes[i].set_title(ax_titles[i], fontsize = 15)
    axes[i].set_xlim(min_max)
    axes[i].set_ylim(min_max)
    axes[i].axhline(y=0, color='k')
    axes[i].axvline(x=0, color='k')

filename = savepath +\
           'Figure2.pdf'
plt.savefig(filename,
            dpi=300)

