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
            dataframe with metabolite data
        QT_pad: clean.QT_pad
            dataframe with imaging data and covariates

        Attributes
        ----------
        data: list pd.DataFrame
            merged dataframe split by sex
        sexes: list
            list of the order of sexes
        phenotypes: list
            phenotype names
        covariates: list
            covariate names
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

        #### Split by sex
        dats = []
        for i in range(2):
            bools = dat['PTGENDER'] == self.sexes[i]
            dats.append(dat.loc[bools])
        
        for i in range(2):
            dats[i].drop('PTGENDER', 
                         axis=1)
        
        self.data = dats

    def run_ewas(self):
        '''
        Run an EWAS for each phenotype

        Returns
        ----------
        results: list of pd.DataFrame
            list of results based on sex by phenotype order 
        '''
        results = []
        for s in range(len(self.sexes)):
            for p in self.phenotypes:
                drop_phenos  = [m for m in self.phenotypes \
                                if m != p]
                data_temp = self.data[s].drop(columns=drop_phenos)
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
            ws = []
            for l in range(len(indices[i])):
                res_temp = self.results[indices[i][l]]
                w = np.array(1 / np.power(\
                            res_temp['SE'], 2))
                ws.append(w)
            meta_se = np.sqrt( 1 / sum(ws) )
            up_term = np.zeros(meta_se.shape)
            for m in range(len(indices[i])):
                temp = np.array(res_temp['Beta'] * ws[m])
                up_term = up_term + temp
                
            meta_beta = up_term / sum(ws)
            zeta = meta_beta / meta_se
            pval = 2 * stats.norm.cdf(-abs(zeta))

            final_dat = pd.DataFrame(meta_se, 
                                    index=res_temp.index, 
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
        
            t2 = np.sqrt(np.power(np.array(res_temp1['SE']), 2) + \
                         np.power(np.array(res_temp2['SE']), 2))
        
            zdiff = t1 / t2
            pval  = 2*stats.norm.cdf(-abs(zdiff))
            pval_bf = multipletests(pval)[1]
            final_dat = pd.DataFrame(zdiff, 
                                     index=res_temp1.index, 
                                     columns=['Z'])
            final_dat['pvalue'] = pval
            final_dat['pvalue_bonferroni'] = pval_bf
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
        indices = [[0,2], [1,3]]
        filter_t = 10 ** -5
        diff_t   = 0.05/55

        for i in range(len(self.results_diff)):
            dat_female = self.results[indices[i][0]]
            dat_male   = self.results[indices[i][1]]
            # Total difference
            total_diff = self.results_diff[i]['pvalue'] < diff_t

            # Filtering first
            overall_filter = self.results_meta[i]['pvalue'] < filter_t
            if sum(overall_filter) > 0:
                bonf_filter_t  = 0.05 / sum(overall_filter)
            else:
                bonf_filter_t  = 0.05
            filter_diff = self.results_diff[i]['pvalue'] < bonf_filter_t

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

        Returns:
        results_diff: pd.DataFrame
            dataframe of the sex difference plus the category 
        '''
        b_alpha = 0.05/55
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
            final_dat['difference_type'] = \
                        self.results_diff[i]['difference_type']

            final_dat['difference_arnold'] = \
                        self.results_diff[i]['difference_arnold']

            name = 'Preliminary_results_' +\
                   self.phenotypes[i] +\
                   '.csv'
            filename = os.path.join(savepath,
                                    name)
            final_dat.to_csv(filename)
 
        