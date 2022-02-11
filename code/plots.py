import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from matplotlib import cm

def female_male_scatter(results,
                        ax,
                        component,
                        modules=None):
    '''
    Plot a scatter with male and female effect sizes in each axis.
    Plot different groups with its corresponding color

    Parameters
    ----------
    results: pd.DataFrame
        result dataframe with 'Beta_female', 'Beta_male', 
        'pvalue_diff', 'Modules', and 'difference_arnold' columns
    ax: plt.axes
        ax to use for matplotlib
    component: str
        component to use, or 'All'
    modules: list of str, str, or None
        module colors to highlight
    '''
    color_pastel = cm.get_cmap('Set2')
    axis_names = ['Effect in females',
                  'Effect in males']
    if component == 'All':
        title_name = 'Single metabolite differences\nacross components'
        dat = results
    else:
        title_name = 'Single metabolite differences in\n' +\
                      component
        component_to_use = results.index.\
                           get_level_values('Outcome') == component
        dat = results[component_to_use]

    sex_different = (dat['difference_arnold'] == 'Sex-specific') |\
                    (dat['difference_arnold'] == 'Heterogeneous')
    everything_else = ~sex_different

    # Divide sex_different by group of metabolites
    meta_names = dat.index.get_level_values('Variable')

    PCs = meta_names.str.match(pat = 'PC.a[ae]') & sex_different
    lysoPCs = meta_names.str.match(pat = 'lysoPC') & sex_different
    acyl = meta_names.str.match(pat = 'C[0-9]') & sex_different
    sphyn = meta_names.str.match(pat = 'SM.') & sex_different
    amino = meta_names.str.match(pat = '[A-Z][a-z]{2}$') & sex_different
    amines = ~(amino + PCs + lysoPCs + acyl + sphyn) & sex_different

    group_names = ['PCs',
                   'lysoPC',
                   'SMs',
                   'Acylcarnitines',
                   'Amino acids',
                   'Biogenic amines']
    
    modules_to_highlight = []
    mod_colors = []
    if modules is not None:
        for mod in modules:
            temp_mod = dat['Modules'] == mod
            sex_different = (sex_different) &\
                            (~temp_mod)
            everything_else = (everything_else) &\
                              (~temp_mod)
            modules_to_highlight.append(temp_mod)
        mod_colors = modules        
    
    groups = [everything_else,
              PCs,
              lysoPCs,
              sphyn,
              acyl,
              amino,
              amines]
    groups.extend(modules_to_highlight)

    colors = []
    for i in range(len(groups)):
        c = color_pastel.reversed()(i)
        colors.append(c)

    colors = colors + mod_colors

    max_value = max(abs(dat.loc[:,['Beta_female',
                                   'Beta_male'] ]).max())
    beta_female = dat['Beta_female']
    beta_male   = dat['Beta_male']
    pvalue_diff = -np.log(dat['pvalue_diff'])*20
    max_value_round = round(max_value + 0.1, 2)
    min_max   = (-max_value_round,
                  max_value_round)

    patches = []
    for i, log in enumerate(groups):
        log = list(log)
        if i == 0:
            al = 0.2
        else:
            al = 0.8
            p = mlines.Line2D([], [], 
                              color=colors[i],
                              marker='o',
                              label=group_names[i-1],
                              ms=10,
                              ls='')
            patches.append(p)
        ax.scatter(beta_female[log],
                   beta_male[log],
                   s=pvalue_diff[log],
                   color=colors[i],
                   alpha=al)
    ax.set_xlabel(axis_names[0], fontsize = 12)
    ax.set_ylabel(axis_names[1], fontsize = 12)
    ax.set_title(title_name, fontsize = 14)
    ax.set_xlim(min_max)
    ax.set_ylim(min_max)
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')
    ax.legend(handles=patches)

def female_male_forest(results,
                       colors,
                       ax):
    '''
    Sex-stratified forest plot of metabolites modules

    Parameters
    ----------
    results: pd.DataFrame
        metabolite module results, with Component 1 on top of the table
    colors: matplotlib.colors.ListedColormap
        colormap to use. Only takes the first to colors for female and male
        in that order
    ax: plt.axes
        ax to use for matplotlib
    '''
    sexes = ['female',
             'male']
    tick_names = results.index.get_level_values('Variable')
    tick_names = [s.strip('ME') for s in tick_names]

    ax.set_title('Differences metabolite modules')
    ax.vlines(0,
              -0.5,
              len(tick_names),
              color='black',
              alpha=0.8,
              linestyles='solid')
    ax.set_yticks(list(range(len(results))))
    ax.set_yticklabels(tick_names)

    # Annotate component names
    component_names = list(results.index.get_level_values('Outcome').unique())
    component_names.reverse()
    step = 0.9/len(component_names)
    pos = step/2
    for comp in component_names:
        ax.annotate(comp,
                    xy=(1.03, pos),
                    xycoords='axes fraction',
                    rotation=90)
        pos = pos + step

    patches = []
    for sex in sexes:
        if sex == 'female':
            sep = 0
            color = colors(0)
        else:
            sep = 0.3
            color = colors(1)
        beta_name = 'Beta_' + sex
        se_name   = 'SE_' + sex
        interval  = st.norm.interval(alpha=0.95,
                                     loc=results[beta_name],
                                     scale=results[se_name])
        for y_row in range(len(interval[0])):
            c = y_row + sep
            ax.scatter(results[beta_name][y_row],
                       c,
                       color=color)
            ax.plot([round(interval[0][y_row],3),
                     round(interval[1][y_row],3)],
                    [c,c],
                     color=color)
        p = mlines.Line2D([], [], 
                          color=color,
                          marker='o',
                          label=sex,
                          ms=10,
                          ls='-')
        patches.append(p)
        
    ax.legend(handles=patches)

def score_plot(ax,
               qtpad,
               components:list=['Component 1', 'Component 2']):
    '''
    Generate a scatter plot from PLS-DA scores

    Parameters
    ----------
    ax: ax
        Matplotlib ax to use
    qtpad: Qtpad
        qtpad class with plsda scores generated
    components: list of str
        Name of Components to plot. Must match score column names
    '''
    percents = _get_percents(qtpad,
                             components)
    color_pastel = cm.get_cmap('Set2')
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
    for c, target in enumerate(targets):
        indicesToKeep = qtpad.data['DX.bl'] == target
        ax.scatter(qtpad.scores.loc[indicesToKeep,
                                    components[0]],
                   qtpad.scores.loc[indicesToKeep,
                                    components[1]],
                   color = color_pastel(c),
                   s = 50,
                   alpha = 0.5)
    ax.legend(targets)

def weight_plot(ax,
                qtpad,
                weights:list=['Weight 1', 'Weight 2'],
                offset:dict=None):
    '''
    Plot the weights of phenotypes and diagnosis

    Parameters
    ----------
    ax: ax
        Matplotlib ax to use
    qtpad: Qtpad
        qtpad class with plsda scores generated
    weights: list of str
        Name of Weights to plot. Must match weight column names
    offset: dict
        Phenotype offset dict to change annotation possition. 
        E.g. {'Ventricles': (-10, 10)}
    '''
    percents = _get_percents(qtpad,
                             weights)
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
        xy_offset = (-50, 4)
        if offset is not None:
            if segment in offset.keys():
                xy_offset = offset[segment]

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

def vip_plot(ax,
             qtpad):
    '''
    Variable importance in projection (VIP) Plot

    Parameters
    ----------
    ax: ax
        Matplotlib ax to use
    qtpad: Qtpad
        Qtpad class with PLS-DA applied
    '''
    font_ax_title = 16
    pos = np.arange(len(qtpad.vips))
    order = qtpad.vips.argsort()[::-1]
    ax.set_title('Variable importance in projection\n(VIP)',
                 fontsize = font_ax_title)
    ax.barh(pos,
            qtpad.vips[order],
            align='center')
    ax.set_yticks(pos)
    ax.set_yticklabels(np.array(qtpad.phenotypes)[order],
                       rotation=60)

def scree_plot(ax,
               qtpad):
    '''
    Scree Plot

    Parameters
    ----------
    ax: ax
        Matplotlib ax to use
    qtpad: Qtpad
        Qtpad class with PLS-DA applied
    '''
    font_ax_title = 16
    pcs = np.arange(len(qtpad.x_variance_explained))
    ax.plot(pcs, qtpad.x_variance_explained)
    ax.set_title('Scree Plot',
                 fontsize = font_ax_title)
    ax.set_ylabel('Variance Explained')
    ax.set_xticks(pcs)
    ax.set_xticklabels(list(qtpad.scores.columns),
                       rotation=60)

def _get_percents(qtpad,
                  components:list):
    '''
    Get the percents of variance explained of the components or weights
    out of component or weight names
    
    Parameters
    ----------
    qtpad: Qtpad
        qtpad class with PLS-DA applied
    components: list of str
        list of column names (components or weights)
    
    Returns
    ----------
    percents: list of str
        list of percents explained
    '''
    percents = []
    for c in components:
        i = [int(s) for s in c.split() if s.isdigit()][0] - 1
        per = '(' + \
                str(round(qtpad.x_variance_explained[i] * 100)) + \
                '%)'
        percents.append(per)
    return(percents)
    