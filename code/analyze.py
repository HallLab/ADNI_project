import os
import clean
import clarite
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

class ADNI:
    '''
    ADNI class that integrates imaging data, metabolites, and covariates
    '''
    def __init__(self,
                 metabolite_type:str='p180',
                 modules:bool=False,
                 phenotypes_pls:bool=False):
        '''
        Initiate the class
        
        Parameters
        ----------
        metabolite_type: str
            type of metabolite data to load, either p180 or nmr
        modules: bool
            whether to load wgcna modules or single metabolites
        phenotypes_pls: bool
            whether to load the pls scores instead of individuals phenotypes

        Attributes
        ----------
        data: pd.DataFrame
            dataframe containing the values
        metabolite_names: list of str
            metabolite names
        phenotype_names: list of str
            phenotype names
        covariate_names: list of str
            covariate names
        modules: bool
            whether it is modules or not
        filename: str
            name of the metabolite_type and module for future saving
        '''
        # Define filename
        if modules:
            name = '_modules'
        else:
            name = ''     
        self.filename = 'results_' + \
                        metabolite_type + \
                        name
        self.modules = modules
        
        respath = '../results/'
        
        # Load phenotypes and covariates
        qtpad = clean.QT_pad()
        covs   = qtpad.data.loc[:,qtpad.covariates]
        self.covariate_names = qtpad.covariates
        if phenotypes_pls:
            phenos = pd.read_csv('../results/pheno_pls_components.csv').\
                        set_index('RID')
            self.phenotype_names = ['Component 1',
                                    'Component 2']
        else:
            phenos = qtpad.data.loc[:,qtpad.phenotypes]
            self.phenotype_names = qtpad.phenotypes
            
        # Load metabolite data
        if metabolite_type == 'p180':
            if modules:
                file_metabolite = respath +\
                                  'eigenmetabolites_p180.csv'
            else:
                file_metabolite = respath +\
                                  'p180_cleaned.csv'
        elif metabolite_type == 'nmr':
            if modules:
                file_metabolite = respath +\
                                  'eigenmetabolites_nmr.csv'
            else:
                file_metabolite = respath +\
                                  'nmr_cleaned.csv'
        else:
            print('Please provide an appropriate metabolite_type (p180 or nmr)')
            return()
        
        metabolites = pd.read_csv(file_metabolite).\
                                  set_index('RID')
        self.metabolite_names = list(metabolites.columns)

        # Three way merge
        one_merge = covs.merge(phenos,
                               left_on='RID',
                               right_on='RID')
        final_merge = one_merge.merge(metabolites,
                                      left_on='RID',
                                      right_on='RID')

        self.data = final_merge         

    def normalize_data(self):
        '''
        Zscore normalize the phenotype, metabolites, and covariate data
        stratified by sex

        Returns
        ----------
        normalized_data: pd.DataFrame
                         normalized data
        '''
        names_to_normalize = self.covariate_names +\
                             self.phenotype_names +\
                             self.metabolite_names
        remove_names = ['PTGENDER',
                        'APOE4']
        for name in remove_names:
            names_to_normalize.remove(name)

        data_to_normalize = self.data[names_to_normalize]
        qtpad = clean.QT_pad()
        normalized_data = clean.zscore_normalize(data_to_normalize,
                                                 qtpad)
        self.data[names_to_normalize] = normalized_data

    def stratified_association(self):
        '''
        Run a stratified association for each phenotype

        Returns
        ----------
        results: list of pd.DataFrame
            list of results based on sex by phenotype order 
        '''
        results = []
        data    = []
        sexes   = ['Female',
                   'Male']
        if self.modules:
            # Don't analyze grey module
            module_names = self.metabolite_names.copy()
            if 'MEgrey' in module_names:
                module_names.remove('MEgrey')
            for s in range(len(sexes)):
                keep_cols = module_names + \
                            self.phenotype_names + \
                            self.covariate_names
                temp_data = self.data.loc[:,keep_cols]
        else:
            temp_data = self.data
        
        for s in sexes:
            bool_sex = temp_data['PTGENDER'] == s
            d = temp_data.loc[bool_sex]
            d.drop('PTGENDER', 
                   axis=1,
                   inplace=True)
            data.append(d)

        categorical_vars = 'APOE4'
    
        for i in range(len(data)):                
            for p in self.phenotype_names:
                drop_phenos  = [m for m in self.phenotype_names \
                                if m != p]
                data_temp = data[i].drop(columns=drop_phenos)
                # Change RID to ID because clarite
                data_temp.index.rename('ID',
                                       inplace=True)
                # Categorize data
                continuous_vars = list(data_temp.columns)
                continuous_vars.remove(categorical_vars)
                clarite.modify.make_categorical(data_temp,
                                                categorical_vars)
                clarite.modify.make_continuous(data_temp,
                                               continuous_vars)
                covs = self.covariate_names.copy()
                covs.remove('PTGENDER')
                r = clarite.analyze.ewas(p,
                                         covs,
                                         data_temp)
                results.append(r)

        for i in range(len(results)):
            results[i] = results[i].sort_index()

        self.results = results
        self.meta_analyze()
        self.sex_diff_test()
        self.categorize_sex_diff()
        self.categorize_arnold20()

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
            final_dat = pd.DataFrame(zdiff)
            final_dat.columns = pd.Index(['Z'])
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
                    savepath:str):
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
            final_dat['z_diff'] = self.results_diff[i]['Z']
            final_dat['difference_type'] = \
                        self.results_diff[i]['difference_type']

            final_dat['difference_arnold'] = \
                        self.results_diff[i]['difference_arnold']

            name = self.filename +\
                   '_' +\
                   self.phenotype_names[i] +\
                   '.csv'
            fullname = os.path.join(savepath,
                                    name)
            final_dat.to_csv(fullname)
 