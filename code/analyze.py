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
        metabolite_dict: pd.DataFrame
            dataframe with extra notes and information of metabolites
            (only if modules is False)
        metabolite_modules: pd.DataFrame
            dataframe with metabolite module assignment
            (only if modules is False)
        metabolite_names: list of str
            metabolite names
        phenotype_names: list of str
            phenotype names
        covariate_names: list of str
            covariate names
        modules: bool
            whether it is metabolite wgcna modules or not
        filename: str
            name of the metabolite_type and module for future saving
        '''
        # Define filename
        respath = '../results/'
        datapath = '../data/'
        if modules:
            name = '_modules'
        else:
            name = ''
            usecols = ['FLDNAME',
                       'TEXT',
                       'NOTES']
            if metabolite_type == 'p180':
                dict_names = [datapath + 
                              'ADMCDUKEP180FIAADNI2GO_DICT.csv',
                              datapath +
                              'ADMCDUKEP180UPLCADNI2GO_DICT.csv']
                metabolite_dict = []
                for n in dict_names:
                    dat = pd.read_csv(n,
                                      usecols=usecols).\
                             set_index('FLDNAME')
                    metabolite_dict.append(dat)
                metabolite_dict = pd.concat(metabolite_dict,
                                            axis=0)
                metabolite_dict['KEGG'] = metabolite_dict['NOTES'].\
                                            str.\
                                            extract(r"(KEGG = C\d+)" )
                metabolite_dict['HMDB'] = metabolite_dict['NOTES'].\
                                            str.\
                                            extract(r"(HMDB = HMDB\d+)" )
                to_replace = ['KEGG',
                              'HMDB']
                for rep in to_replace:
                    metabolite_dict[rep] = metabolite_dict[rep].\
                                            str.\
                                            replace(rep +
                                                    ' = ',
                                                    '')
                metabolite_dict.drop('NOTES',
                                     axis=1,
                                     inplace=True)
                module_colors_name = 'module_colors_p180.csv'
            elif metabolite_type == 'nmr':
                dict_names = datapath +\
                             'ADNINIGHTINGALE2_DICT.csv'
                metabolite_dict = pd.read_csv(dict_names,
                                              usecols=usecols).\
                                     set_index('FLDNAME')
                metabolite_dict.rename(columns={'NOTES': 'Class'},
                                       inplace=True)
                module_colors_name = 'module_colors_nmr.csv'

            rename_cols = {'TEXT': 'Description'}
            metabolite_dict.index.rename('Variable',
                                         inplace=True)
            metabolite_dict.rename(columns=rename_cols,
                                   inplace=True)
            metabolite_dict = metabolite_dict.drop_duplicates()
            self.metabolite_dict = metabolite_dict

            rename_meta_cols = {0: 'Modules',
                                1: 'Variable'}
            metabolite_modules = pd.read_csv(respath + module_colors_name,
                                             header=None).\
                                    rename(columns=rename_meta_cols).\
                                    set_index('Variable')
            self.metabolite_modules = metabolite_modules
           
        self.filename = 'results_' + \
                        metabolite_type + \
                        name
        self.modules = modules
        
        # Load phenotypes and covariates
        qtpad = clean.QT_pad()
        covs  = qtpad.data.loc[:,qtpad.covariates]
        self.covariate_names = qtpad.covariates
        if phenotypes_pls:
            phenos = pd.read_csv('../results/pheno_pls_components.csv').\
                        set_index('RID')
            self.phenotype_names = list(phenos.columns)
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
                r = clarite.analyze.association_study(data = data_temp,
                                                      outcomes = p,
                                                      covariates = covs)
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
        indices = []
        half_length = len(self.results) // 2
        for l in range(half_length):
            ind = [l, l+half_length]
            indices.append(ind)
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
        indices = []
        half_length = len(self.results) // 2
        for l in range(half_length):
            ind = [l, l+half_length]
            indices.append(ind)
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

    def categorize_sex_diff(self):
        '''
        Categorize the results from the sex difference.
        Be sure to run sex_diff_test and meta_analysis before

        Returns
        ----------
        results_diff: pd.DataFrame
            dataframe of the sex difference plus the category 
        '''
        indices = []
        half_length = len(self.results) // 2
        for l in range(half_length):
            ind = [l, l+half_length]
            indices.append(ind)
        filter_t = 10**-5
        if self.modules:
            n_modules = len(self.metabolite_names)-1
            diff_t = 0.05/n_modules
        else:
            m = _estimate_effective_tests(self.data[self.metabolite_names])
            diff_t = 0.05/m
            print('There are ' + 
                  str(m) + 
                  ' independent tests')

        for i in range(len(self.results_diff)):
            dat_female = self.results[indices[i][0]]
            dat_male   = self.results[indices[i][1]]
            # Total difference
            total_diff = self.results_diff[i]['pvalue'] < diff_t

            # Filtering first
            overall_filter = self.results_meta[i]['pvalue'] < filter_t
            if sum(overall_filter) > 0:
                filtered_metabolites = list(overall_filter[\
                                            overall_filter].\
                                                index.\
                                                get_level_values(0))
                filter_m = _estimate_effective_tests(self.data\
                                                    [filtered_metabolites])
                print('There are ' + 
                      str(filter_m) + 
                      ' filtered independent tests')
                bonf_filter_t = 0.05 / filter_m
                filter_diff   = self.results_diff[i].\
                                     loc[overall_filter,
                                         'pvalue'] < bonf_filter_t
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

    def categorize_arnold20(self):
        '''
        Categorize based on criteria from Arnold et al 2020 paper

        Returns
        ----------
        results_diff: pd.DataFrame
            dataframe of the sex difference plus the category 
        '''
        if self.modules:
            n_modules = len(self.metabolite_names)-1
            b_alpha = 0.05/n_modules
        else:
            if 'p180' in self.filename:
                b_alpha   = 0.05/50
            elif 'nmr' in self.filename:
                b_alpha   = 0.05/44

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
        indices = []
        half_length = len(self.results) // 2
        for l in range(half_length):
            ind = [l, l+half_length]
            indices.append(ind)
        final_dat = []
        
        for i in range(len(indices)):
            f = indices[i][0]
            m = indices[i][1]
            dat = self.results[f][['Beta',
                                   'SE',
                                   'pvalue']].copy()
            dat.columns = ['Beta_female',
                           'SE_female',
                           'pvalue_female']
            dat['Beta_male'] = self.results[m]['Beta']
            dat['SE_male']   = self.results[m]['SE']
            dat['pvalue_male'] = self.results[m]['pvalue']

            dat['Beta_total'] = self.results_meta[i]['Beta']
            dat['SE_total']   = self.results_meta[i]['SE']
            dat['pvalue_total'] = self.results_meta[i]['pvalue']

            dat['pvalue_diff'] = self.results_diff[i]['pvalue']
            cols_to_round = ['Beta_female',
                             'Beta_male',
                             'Beta_total',
                             'SE_female',
                             'SE_male',
                             'SE_total']
            cols_to_sc_round = ['pvalue_female',
                                'pvalue_male',
                                'pvalue_total',
                                'pvalue_diff']
            dat[cols_to_round] = dat[cols_to_round].round(3)
            for col in cols_to_sc_round:
                dat[col] = dat[col].\
                            apply(lambda x: '{:.1e}'.format(x))
            dat['difference_type'] = \
                        self.results_diff[i]['difference_type']

            dat['difference_arnold'] = \
                        self.results_diff[i]['difference_arnold']
            if self.modules == False:
                dat = pd.merge(dat,
                               self.metabolite_dict,
                               left_index=True,
                               right_index=True)
                dat = pd.merge(dat,
                               self.metabolite_modules,
                               left_index=True,
                               right_index=True)
            final_dat.append(dat)

        final_dat = pd.concat(final_dat,
                              axis=0)

        name = self.filename +\
               '.csv'
        fullname = os.path.join(savepath,
                                name)
        final_dat.to_csv(fullname)
 
def _estimate_effective_tests(data):
    '''
    Estimate number of effective tests based on 
    Li & Ji, 2005

    Parameters
    ----------
    data: pd.DataFrame
        dataframe with values to estimate effective tests
    '''
    # Get correlation
    r = data.corr()
    # Get eigenvalues from r
    w, v = np.linalg.eig(r)
    w = abs(w)
    # Estimate number of effective tests
    m = sum((w >= 1) + (w - np.floor(w)))
    m = np.floor(m)

    return(int(m))
        