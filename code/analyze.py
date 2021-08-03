import os
import clarite
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

class ADNI_sex:
    '''
    ADNI sex class that integrates imaging data and metabolites 
    '''
    def __init__(self,
                 metabolites,
                 QT_pad):
        '''
        Initiate the class and split by sex

        Parameters
        ----------
        metabolites: pd.DataFrame
            dataframe with metabolite data, with or withouth WGCNA modules
        QT_pad: clean.QT_pad
            dataframe with imaging data and covariates

        Attributes
        ----------
        data: list pd.DataFrame
            merged dataframe split by sex
        sexes: list of str
            list of the order of sexes
        phenotypes: list of str
            phenotype names
        covariates: list of str
            covariate names
        module_colors: list of str
            module colors identified with 'ME'
        '''
        self.sexes = ['Female',
                      'Male']
        self.phenotypes = QT_pad.phenotypes
        self.covariates = QT_pad.covariates

        dat = pd.merge(QT_pad.data, 
                       metabolites, 
                       on='RID')
        dat.index.rename('ID',
                         inplace=True)

        # Identify module names
        cols = dat.columns.to_list()
        self.module_colors = [v for v in cols if v.startswith('ME')]
        
        #### Split by sex
        dats = []
        for i in range(2):
            bools = dat['PTGENDER'] == self.sexes[i]
            dats.append(dat.loc[bools])
        
        for i in range(2):
            dats[i].drop('PTGENDER', 
                         axis=1)
        
        self.data = dats

    def zscore_normalize_eigenmetabolites(self):
        '''
        Zscore normalize the eigenmetabolites stratified by sex

        Returns
        ----------
        data: list of pd.Dataframe
            list of data with normalized eigenmetabolites
        '''
        if not self.module_colors:
            print('There are no eigenmetabolites in data')
        else:
            for i in range(len(self.data)):
                self.data[i].loc[:,self.module_colors] = \
                     self.data[i].loc[:,self.module_colors].\
                         apply(stats.zscore,
                               nan_policy='omit')

    def run_ewas(self,
                 on_modules:bool=False):
        '''
        Run an EWAS for each phenotype

        Parameters
        ----------
        on_modules: bool
            Whether to run the EWAS using single metabolites or eigenmetabolites

        Returns
        ----------
        results: list of pd.DataFrame
            list of results based on sex by phenotype order 
        '''
        results = []
        data    = []
        if on_modules:
            if not self.module_colors:
                print('There are no eigenmetabolites in data')
            else:
                # Don't analyze grey module
                modules = self.module_colors.copy()
                if 'MEgrey' in modules:
                    modules.remove('MEgrey')
                for s in range(len(self.sexes)):
                    keep_cols = modules + \
                                self.phenotypes + \
                                self.covariates
                    d = self.data[s].loc[:,keep_cols]
                    data.append(d)
        else:
            for s in range(len(self.sexes)):
                d = self.data[s].drop(columns=self.module_colors)
                data.append(d)
        
        for i in range(len(data)):                
            for p in self.phenotypes:
                drop_phenos  = [m for m in self.phenotypes \
                                if m != p]
                data_temp = data[i].drop(columns=drop_phenos)
                clarite.modify.categorize(data_temp)
                r = clarite.analyze.ewas(p,
                                         self.covariates,
                                         data_temp)
                results.append(r)

        for i in range(len(results)):
            results[i] = results[i].sort_index()

        self.results = results

    def meta_analyze(self):
        '''
        For each phenotypes meta-analyze the results
        and generate new parameters
    
        Returns
        ----------
        results_meta: list of pd.Dataframe
            list of dataframes with new meta analyzed parameters
        '''
        indices = [[0,2], [1,3]]
        results_meta = []
        
        for i in range(len(indices)):
            res_temp1 = self.results[indices[i][0]]
            res_temp2 = self.results[indices[i][1]]
            meta_se = np.sqrt( 1 / 
                  ( ( 1 / np.power(res_temp1['SE'], 2) ) + \
                    ( 1 / np.power(res_temp2['SE'], 2) ) ) )
            
            meta_beta = ( (res_temp1['Beta'] /
                          np.power(res_temp1['SE'],2)) +
                          (res_temp2['Beta'] /
                          np.power(res_temp2['SE'],2)) ) / \
                        ( (1/np.power(res_temp1['SE'],2)) + 
                          (1/np.power(res_temp2['SE'],2)) ) 

            zeta = meta_beta / meta_se
            pval = 2 * stats.norm.cdf(-abs(zeta))

            final_dat = pd.DataFrame(meta_se, 
                                    index=res_temp1.index, 
                                    columns=['SE'])
            final_dat['Beta']   = meta_beta
            final_dat['pvalue'] = pval
            results_meta.append(final_dat)

        self.results_meta = results_meta
        
    def sex_diff_test(self):
        '''
        Test for sex differences
    
        Returns
        ----------
        results_diff: pd.Dataframe
            dataframe containing the beta, SE and pvalues from
            the two sexes, and their difference test
        '''
        indices = [[0,2], [1,3]]
        results_diff = []

        for i in range(len(indices)):
            res_temp1 = self.results[indices[i][0]]
            res_temp2 = self.results[indices[i][1]]
            t1 = np.array(res_temp1['Beta']) - \
                 np.array(res_temp2['Beta'])
        
            t2 = np.sqrt(np.power(res_temp1['SE'], 2) + \
                         np.power(res_temp2['SE'], 2))
        
            zdiff = t1 / t2
            pval  = 2*stats.norm.cdf(-abs(zdiff))
            final_dat = pd.DataFrame(zdiff, 
                                     index=res_temp1.index, 
                                     columns=['Z'])
            final_dat['pvalue'] = pval
            results_diff.append(final_dat)
        self.results_diff = results_diff

    def categorize_sex_diff(self,
                            on_modules:bool=False):
        '''
        Categorize the results from the sex difference.
        Be sure to run sex_diff_test and meta_analysis before

        Parameters
        ----------
        on_modules: bool
            whether the categorization is done on module results

        Returns
        ----------
        results_diff: pd.DataFrame
            dataframe of the sex difference plus the category 
        '''
        indices = [[0,2], [1,3]]
        filter_t = 0.05
        if on_modules:
            n_modules = len(self.module_colors)-1
            diff_t = 0.05/n_modules
        else:
            diff_t   = 0.05/55

        for i in range(len(self.results_diff)):
            dat_female = self.results[indices[i][0]]
            dat_male   = self.results[indices[i][1]]
            # Total difference
            total_diff = self.results_diff[i]['pvalue'] < diff_t

            # Filtering first
            overall_filter = self.results_meta[i]['pvalue'] < filter_t
            if sum(overall_filter) > 0:
                bonf_filter_t = 0.05 / sum(overall_filter)
                filter_diff   = self.results_diff[i]['pvalue'] < bonf_filter_t
            else:
                filter_diff = self.results_diff[i]['pvalue'] > 100 # All false

            significants = total_diff | \
                           filter_diff

            #### CLASSIFICATION
            # 1. Significant, and both female and male 
            # are significant, and betas opposite
            both_sign  = (dat_female['pvalue'] < 0.05) & \
                         (dat_male['pvalue'] < 0.05)
            opposite_direction = dat_female['Beta'] * \
                                 dat_male['Beta'] < 0

            keep_qual = significants & \
                        both_sign & \
                        opposite_direction

            # 2. Overall nominal significance, zdiff significance bonferroni, both significant and same direction
            same_direction   = dat_female['Beta'] * \
                               dat_male['Beta'] > 0
            keep_quant = significants & \
                         same_direction & \
                         both_sign

            # 3. Overall nominal significance, zdiff significance boferroni, only one significant
            one_sig  = ((dat_female['pvalue'] < 0.05 ) & \
                        (dat_male['pvalue'] > 0.05 ) ) | \
                       ((dat_female['pvalue'] > 0.05 ) & \
                        (dat_male['pvalue'] < 0.05 ) )
            keep_pure = significants & \
                        one_sig

            # Adding classification
            self.results_diff[i]['difference_type'] = 'None'
            self.results_diff[i].loc[keep_qual,'difference_type'] =\
                'Qualitative'
            self.results_diff[i].loc[keep_quant,'difference_type'] =\
                'Quantitative'
            self.results_diff[i].loc[keep_pure,'difference_type']  =\
                'Pure'

    def categorize_arnold20(self, 
                            on_modules:bool=False):
        '''
        Categorize based on criteria from Arnold et al 2020 paper

        Parameters
        ----------
        on_modules: bool
            whether the categorization is done on module results

        Returns
        ----------
        results_diff: pd.DataFrame
            dataframe of the sex difference plus the category 
        '''
        if on_modules:
            n_modules = len(self.module_colors)-1
            b_alpha = 0.05/n_modules
        else:
            b_alpha   = 0.05/55

        for i in range(len(self.results_diff)):
            f = i + 2 # index for male sex specific result
            filter_1 = self.results_meta[i]['pvalue'] < b_alpha

            bool_1 = self.results[i]['pvalue'] < b_alpha
            bool_2 = self.results[f]['pvalue'] < b_alpha
            filter_2 = bool_1 | bool_2

            bool_3 = self.results[i]['pvalue'] < 0.05
            bool_4 = self.results[f]['pvalue'] < 0.05
            bool_5 = self.results_diff[i]['pvalue'] < 0.05
            filter_3 = (bool_3 | bool_4) & bool_5

            final_filter = filter_1 | filter_2 | filter_3

            self.results_diff[i]['difference_arnold'] = 'None'
            self.results_diff[i].loc[final_filter,
                                     'difference_arnold'] = 'Homogeneous'
            bool_heterogeneous = final_filter & bool_5
            self.results_diff[i].loc[bool_heterogeneous,
                                     'difference_arnold'] = 'Heterogeneous'
            bool_sex_specific = filter_2 & bool_5
            self.results_diff[i].loc[bool_sex_specific,
                                     'difference_arnold'] = 'Sex-specific'
    
    def save_to_csv(self,
                    savepath:str,
                    modules:bool=False):
        '''
        Save file to csv 

        Parameters
        ----------
        savepath: str
            folder path to save the file
        '''
        indices = [[0,2], [1,3]]
        
        for i in range(len(indices)):
            f = indices[i][0]
            m = indices[i][1]
            final_dat = self.results[f][['Beta',
                                         'SE',
                                         'pvalue']].copy()
            final_dat.columns = ['Beta_female', 
                                 'SE_female', 
                                 'pvalue_female']
            final_dat['Beta_male'] = self.results[m]['Beta']
            final_dat['SE_male']   = self.results[m]['SE']
            final_dat['pvalue_male'] = self.results[m]['pvalue']

            final_dat['Beta_total'] = self.results_meta[i]['Beta']
            final_dat['SE_total']   = self.results_meta[i]['SE']
            final_dat['pvalue_total'] = self.results_meta[i]['pvalue']

            final_dat['pvalue_diff'] = self.results_diff[i]['pvalue']
            final_dat['difference_type'] = \
                        self.results_diff[i]['difference_type']

            final_dat['difference_arnold'] = \
                        self.results_diff[i]['difference_arnold']

            if modules:
                name = 'Preliminary_results_' +\
                        self.phenotypes[i] +\
                        '_modules.csv'
            else:
                name = 'Preliminary_results_' +\
                        self.phenotypes[i] +\
                        '.csv'
            filename = os.path.join(savepath,
                                    name)
            final_dat.to_csv(filename)
 
        