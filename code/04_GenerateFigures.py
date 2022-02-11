#### LIBRARIES ####
import clean
import plots
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
                          'CytoscapeInput-edges-brown-blue-yellow.png')
nmr_network  = plt.imread(savepath + 
                          'CytoscapeInput-edges-brown-turquoise.png')


#### Figure 1 : PLS-DA  ####
fig = plt.figure(figsize = (16,16))
fig.tight_layout()
fig.suptitle('PLS-DA', fontsize=25)
space = 0.3
gs = gridspec.GridSpec(3, 3, wspace=space, hspace=space)

font_ax_title = 16
font_axis = 12
axs = []
for row,col in enumerate([[0,1],[0,1],[0,1]]):
    axs.append(fig.add_subplot(gs[row,col[0]]))
    axs.append(fig.add_subplot(gs[row,col[1]]))

axs.append(fig.add_subplot(gs[0:2,2]))
axs.append(fig.add_subplot(gs[2,2]))

## Score Plots ##
components = [['Component 1', 'Component 2'],
              ['Component 3', 'Component 4'],
              ['Component 4', 'Component 5']]

for i,a in enumerate(range(0,6,2)):
    plots.score_plot(axs[a],
                     qtpad,
                     components[i])

## Weight Plots ##
weights = [['Weight 1', 'Weight 2'],
           ['Weight 3', 'Weight 4'],
           ['Weight 4', 'Weight 5']]

offset_dict = [{'Ventricles': (-5, 10),
                'Entorhinal': (-55, -1),
                'Hippocampus': (-75, -5),
                'WholeBrain': (-65, -5)},
               {'Ventricles': (-5, 10),
                'Entorhinal': (-5, 10),
                'MidTemp': (10, 0),
                'Hippocampus': (-75, -5)},
               {'Ventricles': (-60, -5),
                'WholeBrain': (-5, 10),
                'Entorhinal': (-5, 10)}]

for i,a in enumerate(range(1,6,2)):
    plots.weight_plot(axs[a],
                      qtpad,
                      weights[i],
                      offset_dict[i])

## VIP Plot ##
plots.vip_plot(axs[6],
               qtpad)

## Scree Plot ##
plots.scree_plot(axs[7],
                 qtpad)

## Annotate labels ##
axes = [axs[0],
        axs[1],
        axs[6],
        axs[7]]
for i, label in enumerate(('A', 'B', 'C', 'D')):
    axes[i].text(-0.05, 1.1, 
                 label,
                 transform=axes[i].transAxes,
                 fontsize=16,
                 fontweight='bold',
                 va='top',
                 ha='right')

## Save Figure ##
for ext in file_extensions:
    filename = savepath +\
               'Figure1' +\
               ext
    plt.savefig(filename,
                dpi=300)


#### FIGURE 3: P180 platform ####
## Set Figure
fig = plt.figure(figsize = (12,12))
fig.suptitle('Sex differences in metabolite modules\n(P180 platform)',
             fontsize=18)
gs = gridspec.GridSpec(2, 2, wspace=space)
axs = []
axs.append(fig.add_subplot(gs[0:,0]))
for row in range(2):
    axs.append(fig.add_subplot(gs[row,1]))

## Forest Plot ##
plots.female_male_forest(results_p180_modules,
                         color_pastel,
                         axs[0])

## Network image ##
axs[1].imshow(p180_network)
axs[1].axis('off')
axs[1].set_title('Blue, brown, and yellow\nmetabolite modules')

## Scatter Plot Single Metabolites ##
plots.female_male_scatter(results_p180,
                          axs[2],
                          'All')

for i, label in enumerate(('A', 'B', 'C')):
    axs[i].text(-0.05, 1.1, 
                label,
                transform=axs[i].transAxes,
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
## Set Figure
fig = plt.figure(figsize = (12,12))
fig.suptitle('Sex differences in metabolite modules\n(NMR platform)',
             fontsize=18)
gs = gridspec.GridSpec(2, 2, wspace=space)
axs = []
axs.append(fig.add_subplot(gs[0:,0]))
for row in range(2):
    axs.append(fig.add_subplot(gs[row,1]))

## Forest Plot ##
plots.female_male_forest(results_nmr_modules,
                         color_pastel,
                         axs[0])

## Network image ##
axs[1].imshow(nmr_network)
axs[1].axis('off')
axs[1].set_title('Brown and turquoise\nmetabolite modules')

## Scatter Plot Single Metabolites ##
plots.female_male_scatter(results_nmr,
                          axs[2],
                          'All')

for i, label in enumerate(('A', 'B', 'C')):
    axs[i].text(-0.05, 1.1, 
                label,
                transform=axs[i].transAxes,
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
