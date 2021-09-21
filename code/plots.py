import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt

def female_male_scatter(results,
                        modules,
                        component,
                        ax):
    '''
    Plot a scatter with male and female effect sizes in each axis.
    Plot different groups with its corresponding color

    Parameters
    ----------
    results: pd.DataFrame
        result dataframe with 'Beta_female', 'Beta_male', 
        'pvalue_diff', 'Modules', and 'difference_arnold' columns
    modules: list of str
        module colors to highlight
    component: str
        component to use, either Component 1, or Component 2
    ax: plt.axes
        ax to use for matplotlib
    '''
    axis_names = ['Effect in females',
                  'Effect in males']
    title_name = 'Single metabolite differences in ' +\
                 component
    component_to_use = results.index.\
                       get_level_values('Outcome') == component
    dat = results[component_to_use]

    sex_different = (dat['difference_arnold'] == 'Sex-specific') |\
                    (dat['difference_arnold'] == 'Heterogeneous')
    everything_else = ~sex_different

    modules_to_highlight = []
    for mod in modules:
        temp_mod = dat['Modules'] == mod
        sex_different = (sex_different) &\
                        (~temp_mod)
        everything_else = (everything_else) &\
                          (~temp_mod)
        modules_to_highlight.append(temp_mod)     
    
    groups = [everything_else,
              sex_different]
    groups.extend(modules_to_highlight)

    colors = ['gray'] +\
             ['purple'] +\
             modules

    max_value = max(abs(dat.loc[:,['Beta_female',
                                   'Beta_male'] ]).max())
    beta_female = dat['Beta_female']
    beta_male   = dat['Beta_male']
    pvalue_diff = -np.log(dat['pvalue_diff'])*20
    max_value_round = round(max_value + 0.1, 2)
    min_max   = (-max_value_round,
                  max_value_round)

    for i, log in enumerate(groups):
        if colors[i] == 'gray':
            al = 0.4
        else:
            al = 0.7
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
    ax.annotate('Component 2',
                xy=(1.03, .6),
                xycoords='axes fraction',
                rotation=90)
    ax.annotate('Component 1',
                xy=(1.03, .15),
                xycoords='axes fraction',
                rotation=90)

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
    ax.legend(sexes)
    leg = ax.get_legend()
    leg.legendHandles[0].set_color(colors(0))
    leg.legendHandles[1].set_color(colors(1))
