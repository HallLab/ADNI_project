#### LIBRARIES ####
import enum
import clean
import plots
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib import cm

#### PATHS AND DATA ####
savepath = '../results/plots/'
respath  = '../results/'
file_extensions = ['.pdf',
                   '.jpg']

qtpad = clean.QT_pad()
qtpad.PLS_DA()
qtpad.remove_multivariate_outliers()
color_pastel = cm.get_cmap('Set2')
results_p180_modules = pd.read_csv(respath + 'results_p180_modules.csv').\
                          set_index(['Variable',
                                     'Outcome'])
results_nmr_modules = pd.read_csv(respath + 'results_nmr_modules.csv').\
                         set_index(['Variable',
                                    'Outcome'])
results_p180 = pd.read_csv(respath + 'results_p180.csv').\
                  set_index(['Variable',
                             'Outcome'])
results_nmr = pd.read_csv(respath + 'results_nmr.csv').\
                 set_index(['Variable',
                            'Outcome'])
p180_network = plt.imread(savepath + 
                          'CytoscapeInput-edges-red.txt.png')
nmr_network  = plt.imread(savepath + 
                          'CytoscapeInput-edges-green-turquoise.txt.png')


#### Figure 1
## PLSDA 
comp1 = qtpad.scores['Component 1']
extreme_neg = comp1.sort_values().index[:20]
extreme_pos = comp1.sort_values(ascending=False).index[:20]
round(qtpad.data.loc[extreme_neg].mean(),3)
round(qtpad.data.loc[extreme_pos].mean(),3)

comp2 = qtpad.scores['Component 2']
extreme_neg = comp2.sort_values().index[:20]
extreme_pos = comp2.sort_values(ascending=False).index[:20]
round(qtpad.data.loc[extreme_neg].mean(),3)
round(qtpad.data.loc[extreme_pos].mean(),3)

fig = plt.figure(figsize = (12,12))
fig.tight_layout(h_pad=10)
fig.suptitle('PLS-DA', fontsize=18)

font_ax_title = 16
font_axis = 12
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)
percent_1 = '(' + \
            str(round(qtpad.x_variance_explained[0] * 100)) + \
            '%)'
percent_2 = '(' + \
            str(round(qtpad.x_variance_explained[1] * 100)) + \
            '%)'

#### AX 1 ####
ax1.set_xlabel('Component 1 ' + percent_1, 
               fontsize = font_axis)
ax1.set_ylabel('Component 2 ' + percent_2,
               fontsize = font_axis)
ax1.set_title('Scores',
              fontsize = font_ax_title)
targets = ['CN',
           'MCI',
           'AD']
c = 0
for target in targets:
    indicesToKeep = qtpad.data['DX.bl'] == target
    ax1.scatter(qtpad.scores.loc[indicesToKeep,
                                 'Component 1'],
                qtpad.scores.loc[indicesToKeep,
                                 'Component 2'],
                color = color_pastel(c),
                s = 50,
                alpha = 0.6)
    c = c + 1
ax1.legend(targets)

#### AX 2 ####
ax2.set_xlabel('Component 1 ' + percent_1, 
               fontsize = font_axis)
ax2.set_ylabel('Component 2 ' + percent_2,
               fontsize = font_axis)
ax2.set_title('Weights',
              fontsize = font_ax_title)
ax2.scatter(qtpad.x_weights['Weight 1'],
            qtpad.x_weights['Weight 2'])
for segment in qtpad.x_weights.index:
    if segment == 'Ventricles':
        xy_offset = (-10, 10)
    elif segment == 'Entorhinal':
        xy_offset = (-55, -1)
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

#### AX 3 ####
pcs = np.arange(len(qtpad.x_variance_explained))
ax3.plot(pcs, qtpad.x_variance_explained)
ax3.set_title('Scree Plot')
ax3.set_ylabel('Variance Explained')
ax3.set_xticks(pcs)
ax3.set_xticklabels(list(qtpad.scores.columns),
                    rotation=60)

#### AX 4 ####
pos = np.arange(len(qtpad.vips))
order = qtpad.vips.argsort()[::-1]
ax4.set_title('Variable importance in projection (VIP)',
              fontsize = font_ax_title)
ax4.bar(pos,
        qtpad.vips[order],
        align='center',
        width=0.6)
ax4.set_xticks(pos)
ax4.set_xticklabels(np.array(qtpad.phenotypes)[order],
                    rotation=60)

# Annotate labels
axes = [ax1,
        ax2,
        ax3,
        ax4]
for i, label in enumerate(('A', 'B', 'C', 'D')):
    axes[i].text(-0.05, 1.1, 
                 label,
                 transform=axes[i].transAxes,
                 fontsize=16,
                 fontweight='bold',
                 va='top',
                 ha='right')


for ext in file_extensions:
    filename = savepath +\
               'Figure1' +\
               ext
    plt.savefig(filename,
                dpi=300)

#### FIGURE 3: P180 platform ####
# SET FIGURE
fig = plt.figure(figsize = (12,12))
fig.suptitle('Sex differences in P180 platform', fontsize=18)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# AX1: forest plot
plots.female_male_forest(results_p180_modules,
                         color_pastel,
                         ax1)

# AX2: network
ax2.imshow(p180_network)
ax2.axis('off')
ax2.set_title('Red metabolite module')

# AX3: scatter plot Component 1
plots.female_male_scatter(results_p180,
                          ['red'],
                          'Component 1',
                          ax3)
# AX4: scatter plot Component 2
plots.female_male_scatter(results_p180,
                          ['red'],
                          'Component 2',
                          ax4)
# Annotate labels
axes = [ax1,
        ax2,
        ax3,
        ax4]
for i, label in enumerate(('A', 'B', 'C', 'D')):
    axes[i].text(-0.05, 1.1, 
                 label,
                 transform=axes[i].transAxes,
                 fontsize=16,
                 fontweight='bold',
                 va='top',
                 ha='right')

for ext in file_extensions:
    filename = savepath +\
               'Figure3' +\
               ext
    plt.savefig(filename,
                dpi=300)

#### FIGURE 4: NMR platform ####
fig = plt.figure(figsize = (12,12))
fig.suptitle('Sex differences in NMR platform', fontsize=18)
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# AX1: forest plot
plots.female_male_forest(results_nmr_modules,
                         color_pastel,
                         ax1)

# AX2: network
ax2.imshow(nmr_network,
           aspect='auto')
ax2.axis('off')
ax2.set_title('Green and turquoise metabolite modules')

# AX3: scatter plot Component 1
plots.female_male_scatter(results_nmr,
                          ['green', 'turquoise'],
                          'Component 1',
                          ax3)

# AX4: scatter plot Component 2
plots.female_male_scatter(results_nmr,
                          ['green', 'turquoise'],
                          'Component 2',
                          ax4)
# Annotate labels
axes = [ax1,
        ax2,
        ax3,
        ax4]
for i, label in enumerate(('A', 'B', 'C', 'D')):
    axes[i].text(-0.05, 1.1, 
                 label,
                 transform=axes[i].transAxes,
                 fontsize=16,
                 fontweight='bold',
                 va='top',
                 ha='right')
                 
for ext in file_extensions:
    filename = savepath +\
               'Figure4' +\
               ext
    plt.savefig(filename,
                dpi=300)
            
#### SUPP FIG 2 ####
# PLS-DA scatter plots #

def score_plot(ax,
               percents:list,
               components:list=['Component 1', 'Component 2']):
    '''
    Generate a scatter plot from PLS-DA scores

    Parameters
    ----------
    components: list of str
        Name of Components to plot. Must match score column names
    percents: list of str
        Percentages to use for axis names
    ax: ax
        Matplotlib ax to use
    '''
    font_ax_title = 16
    font_axis = 12
    ax.set_xlabel(components[0] + ' ' + percents[0], 
                  fontsize = font_axis)
    ax.set_ylabel(components[1] + ' ' + percents[1],
                  fontsize = font_axis)
    ax.set_title('Scores',
                 fontsize = font_ax_title)
    targets = ['CN',
               'MCI',
               'AD']
    c = 0
    for target in targets:
        indicesToKeep = qtpad.data['DX.bl'] == target
        ax.scatter(qtpad.scores.loc[indicesToKeep,
                                    components[0]],
                   qtpad.scores.loc[indicesToKeep,
                                    components[1]],
                   color = color_pastel(c),
                   s = 50,
                   alpha = 0.6)
        c = c + 1
    ax.legend(targets)

def weight_plot(ax,
                percents:list,
                weights:list=['Weight 1', 'Weight 2'],
                offset:dict=None):
    '''
    Plot the weights of phenotypes and diagnosis

    Parameters
    ----------
    components: list of str
        Name of Components to plot. Must match score column names
    percents: list of str
        Percentages to use for axis names
    ax: ax
        Matplotlib ax to use
    offset: dict
        Phenotype offset dict to change annotation possition. 
        E.g. {'Ventricles': (-10, 10)}
    '''
    font_ax_title = 16
    font_axis = 12
    components = [str(l) for i in weights for l in i.split() if l.isdigit()]
    ax.set_xlabel('Component' + components[0] + ' ' + percents[0],
                  fontsize = font_axis)
    ax.set_ylabel('Component' + components[1] + ' ' + percents[1],
                  fontsize = font_axis)
    ax.set_title('Weights',
                 fontsize = font_ax_title)
    ax.scatter(qtpad.x_weights[weights[0]],
               qtpad.x_weights[weights[1]])
    for segment in qtpad.x_weights.index:
        if offset is not None:
            if segment in offset.keys():
                xy_offset = offset[segment]
        else:
            xy_offset = (-55, 4)
        ax.annotate(segment, 
                    xy=(qtpad.x_weights.loc[segment,
                                           weights[0]],
                        qtpad.x_weights.loc[segment,
                                           weights[1]]),
                    xytext=xy_offset,
                    textcoords='offset points')
    ax.scatter(qtpad.y_weights[weights[0]],
               qtpad.y_weights[weights[1]])
    for group in qtpad.y_weights.index:
        xy_offset = (-10, 5)
        ax.annotate(group, 
                    xy=(qtpad.y_weights.loc[group,
                                           weights[0]],
                        qtpad.y_weights.loc[group,
                                           weights[1]]),
                    xytext=xy_offset,
                    textcoords='offset points')
    

fig = plt.figure(figsize = (12,12))
fig.tight_layout(h_pad=10)
fig.suptitle('PLS-DA Scores', fontsize=18)

ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)
percents = []
for i in range(2, len(qtpad.x_variance_explained)):
    per = '(' + \
            str(round(qtpad.x_variance_explained[i] * 100)) + \
            '%)'
    percents.append(per)

#### AX 1 ####
score_plot(ax1,
           percents[:2],
           ['Component 3', 'Component 4'])

#### AX 2 ####
weight_plot(ax2,
            percents[:2],
            ['Weight 3', 'Weight 4'],
            {'Ventricles': (-10, 10),
            'Entorhinal': (-40, 10)})

#### AX 3 ####
score_plot(ax3,
           percents[1:],
           ['Component 4', 'Component 5'])

#### AX 4 ####
weight_plot(ax4,
            percents[1:],
            ['Weight 4', 'Weight 5'])

for ext in file_extensions:
    filename = savepath +\
               'FigureS2' +\
               ext
    plt.savefig(filename,
                dpi=300)
