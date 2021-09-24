import numpy as np
import pandas as pd
import pingouin as pg
import statsmodels.api as sm

from warnings import simplefilter
from sklearn import linear_model
from sklearn.preprocessing import OneHotEncoder
from sklearn.cross_decomposition import PLSRegression

import scipy.stats as stats

class Metabolites:
    '''
    Raw metabolomics platform class
    '''
    def __init__(self,
                 platform:str):
        '''
        Initiate the class by loading the corresponding metabolomics platform

        Parameters
        ----------
        platform: str
            type of metabolomics platform to load (p180 or nmr)

        Attributes
        ----------
        data: list of pd.DataFrame or pd.DataFrame
            dataframe or list of dataframe with concentration values
        platform: str
            type of metabolomics platform loaded (p180 or nmr)
        fasting: pd.DataFrame
            fasting information

        Attributes (p180 exclusive)
        ----------
        cohort: list of str
            List of the ADNI cohort (ADNI1 or ADNI2GO)
        type: list of str
            List of the type of analysis (UPLC or FIA)
        pool: list of pd.DataFrame
            List of data from sample pools
        data_meta: list of pd.DataFrame
            List of metadata from the participants data
        pool_meta: list of pd.DataFrame
            List of metadata from the pool samples

        Attributes (nmr)
        ----------
        qc_tags: pd.DataFrame
            quality control tags
        '''
        if platform!='p180' and platform!='nmr':
            print('Select an appropriate metabolomics platform (p180 or nmr)')
        else:
            self.platform = platform

        print('-----Loading the ' +
              self.platform + 
              ' metabolomics platform-----')

        #### SETTING PATHS AND FILES ####
        data_path = '../data/'

        #### METABOLOMICS DATA ####
        if self.platform == 'p180':
            file_names = ['ADMCDUKEP180UPLC_01_15_16.csv',
                          'ADMCDUKEP180FIA_01_15_16.csv',
                          'ADMCDUKEP180UPLCADNI2GO.csv',
                          'ADMCDUKEP180FIAADNI2GO.csv']

            self.cohort = ['ADNI1',
                           'ADNI1',
                           'ADNI2GO',
                           'ADNI2GO']
      
            self.type = ['UPLC',
                         'FIA',
                         'UPLC',
                         'FIA']
            
            na_values = ['< LOD',
                         'No Interception',
                         '>Highest CS']
            
            self.pool = []
            self.data = []
            self.lod  = []
            self.data_meta = []
            self.pool_meta = []
            c = 0
            for i in file_names:
                dat = pd.read_csv(data_path + i,
                                  na_values=na_values).\
                         set_index('RID')
                # Carnosine is misspelled in ADNI2GO UPLC
                if i == 'ADMCDUKEP180UPLCADNI2GO.csv':
                    dat = dat.rename(columns={'canosine':'Carnosine'})
                #Divide between pool and proper samples
                dat_pool = dat.loc[999999]
                dat_data = dat.loc[0:99999]
                #Divide between metadata and data (metabolites)
                if self.cohort[c] == 'ADNI1':
                    col_index = list(range(0,7))       
                elif self.cohort[c] == 'ADNI2GO':
                    col_index = list(range(0,24))
                col_index.append(-1)
                self.pool_meta.append(dat_pool.iloc[:,col_index])
                self.data_meta.append(dat_data.iloc[:,col_index])
                dat_data.drop(dat_data.columns[col_index],
                              axis = 1,
                              inplace = True)
                dat_pool.drop(dat_pool.columns[col_index],
                              axis = 1,
                              inplace = True)
                self.data.append(dat_data)
                self.pool.append(dat_pool)
                c = c + 1

            #### LOD VALUES ####
            for i in range(len(self.data)):
                filename = 'P180' + \
                            str(self.type[i]) + \
                            'LODvalues_' + \
                            str(self.cohort[i]) + \
                            '.csv'
                # In lod value ADNI2GO FIA, the bar code plate
                # needs fixing
                if i == 3:
                    dat = pd.read_csv(data_path + filename,
                                      encoding='latin_1')
                    barcode = dat['Plate Bar Code']
                    barcode = barcode.str.split(' ',
                                                expand = True)[2].\
                                      str.replace(pat='/',
                                                  repl='-')
                    dat['Plate Bar Code'] = barcode
                else:
                    dat = pd.read_csv(data_path + filename)
                # Metabolite names in lod don't match those in
                # the data, replace '-', ':', '(', ')' and ' ' with '.'
                old_columns = dat.columns
                new_columns = old_columns.str.replace(pat='-|:|\(|\)| ',
                                                      repl='.')
                dat.columns = new_columns
                # Change metabolite name from Met.So to Met.SO
                if self.type[i] == 'UPLC':
                    dat.rename(columns={'Met.SO':'Met.So'},
                               inplace=True)
                self.lod.append(dat)
    
        else:
            file_names = 'ADNINIGHTINGALE2.csv'
            dat = pd.read_csv(data_path + file_names,
                              na_values='TAG').\
                     set_index(['RID','VISCODE2'])
            dat.drop(['VISCODE',
                      'EXAMDATE',
                      'SAMPLEID',
                      'GLOBAL.SPEC.ID',
                      'update_stamp'],
                     axis=1,
                     inplace=True)
            qc_tag_names = ['EDTA_PLASMA',
                            'CITRATE_PLASMA',
                            'LOW_ETHANOL',
                            'MEDIUM_ETHANOL',
                            'HIGH_ETHANOL',
                            'ISOPROPYL_ALCOHOL',
                            'N_METHYL_2_PYRROLIDONE',
                            'POLYSACCHARIDES',
                            'AMINOCAPROIC_ACID',
                            'LOW_GLUCOSE',
                            'HIGH_LACTATE',
                            'HIGH_PYRUVATE',
                            'LOW_GLUTAMINE_OR_HIGH_GLUTAMATE',
                            'GLUCONOLACTONE',
                            'LOW_PROTEIN',
                            'UNEXPECTED_AMINO_ACID_SIGNALS',
                            'UNIDENTIFIED_MACROMOLECULES',
                            'UNIDENTIFIED_SMALL_MOLECULE_A',
                            'UNIDENTIFIED_SMALL_MOLECULE_B',
                            'UNIDENTIFIED_SMALL_MOLECULE_C',
                            'BELOW_LIMIT_OF_QUANTIFICATION']
            qc_tags = dat.loc[:,qc_tag_names]
            # Remove qc tags columns that have only 0
            remove_qc_cols = qc_tags.columns[qc_tags.sum() == 0]
            qc_tags.drop(remove_qc_cols,
                         axis=1,
                         inplace=True)
    
            dat.drop(qc_tag_names,
                     axis=1,
                     inplace=True)
    
            dat = _keep_baseline(dat)
            dat = dat.droplevel('VISCODE2')
            qc_tags = _keep_baseline(qc_tags)
            qc_tags = qc_tags.droplevel('VISCODE2')

            self.data = dat
            self.qc_tags = qc_tags
        
        #### FASTING DATA ####
        fast_dat = pd.read_csv(data_path + 'BIOMARK.csv', 
                               index_col='RID',
                               na_values=-4)
        #Keep only information from baseline
        fast_dat = fast_dat[fast_dat['VISCODE2'] == 'bl']
        fast_dat = fast_dat['BIFAST']
        #If duplicates, keep the largest observed value
        duplicates_ID = fast_dat.index[\
                        fast_dat.index.duplicated()].unique()
        for i in duplicates_ID:
            val = fast_dat.loc[i].max()
            fast_dat.drop(i,
                          axis='index',
                          inplace=True)
            fast_dat.loc[i] = val
        self.fasting = fast_dat.sort_index()

        #Printing information
        if self.platform == 'nmr':
            n_participants = self.data.shape[0]
            n_metabolites  = self.data.shape[1]
        elif self.platform == 'p180':
            n_participants = self.data[0].shape[0] + \
                             self.data[3].shape[0]
            n_metabolites  = self.data[0].shape[1] + \
                             self.data[1].shape[1]
        print('There are ' + 
              str(n_participants) +
              ' participants, and ' +
              str(n_metabolites) +
              ' metabolites\n')

    def remove_missing_metabolites(self,
                                   cutoff: float = 0.2,
                                   by_cohort: bool = True):
        '''
        Remove metabolites due to missing data greater than cutoff
        
        Parameters
        ----------
        cutoff: float
            Missing data removal cutoff
        by_cohort: bool
            Whether to analyze by cohort, or by merging metabolites
            across cohorts (only applicable in p180 platform)

        Returns
        ----------
        data: pd.DataFrame
            Dataframe with metabolites removed due to missingness
        pool (p180): pd.DataFrame
            Dataframe with metabolites removed due to missingness
        '''
        print('-----Removing metabolites with missing data greater than ' +
              str(cutoff) + '-----')
        if self.platform == 'p180':
            metabolite_list = []
            if by_cohort:
                indices = range(len(self.data))
            else:
                indices = [[0,2], [1,3]]
            for i in indices:
                if by_cohort:
                    data = self.data[i]
                else:
                    data = pd.concat([self.data[i[0]],
                                      self.data[i[1]]])
                remove_met_table = _estimate_delete_missingness(data)
                self._print_metabolites_removed(remove_met_table, i)
                #Remove metabolites from data and pool
                self._remove_metabolites(remove_met_table, i)
                metabolite_list.extend(list(remove_met_table.index))
            metabolite_list = set(metabolite_list)
            n_metabolites   = len(metabolite_list)
            print('There were ' + 
                  str(n_metabolites) + 
                  ' metabolites removed in the p180 platform across cohorts\n')
        elif self.platform == 'nmr':
            remove_met_table = _estimate_delete_missingness(self.data)
            self._print_metabolites_removed(remove_met_table)
            self._remove_metabolites(remove_met_table)
            print('')

    def remove_metabolites_cv(self,
                              cutoff: float = 0.2,
                              by_cohort: bool = True):
        '''
        Compute the coefficient of variation among duplicates or triplicates
        for each metabolite and remove metabolites with CV higher than cutoff.
        Can only be applied to p180 platform

        Parameters
        ----------
        cutoff: float
            CV metabolite removal cutoff.
        by_cohort: bool
            Whether to analyze by cohort, or by merging metabolites
            across cohorts (only applicable in p180 platform)

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with metabolites removed due to high CV.
        pool: pd.Dataframe
            Dataframe with metabolites removed due to high CV.
        '''
        if self.platform == 'nmr':
            print('Cannot compute the CV in the NMR platform')
        else:
            print('-----Removing metabolites with CV values greater than ' +
                  str(cutoff) + '-----')
            metabolite_list = []
            if by_cohort:
                indices = range(len(self.data))
            else:
                indices = [[0,2], [1,3]]
            for i in indices:
                cv_interplate = []
                if by_cohort:
                    data = self.data[i]
                else:
                    data = pd.concat([self.data[i[0]],
                                      self.data[i[1]]])
                duplicates_ID = data.index[\
                                data.index.duplicated()].unique()
                for j in range(len(duplicates_ID)):
                    duplicates_dat = data.loc[duplicates_ID[j],:]
                    cv = duplicates_dat.std() / duplicates_dat.mean()
                    cv_interplate.append(cv)
                cv_interplate = pd.DataFrame(pd.DataFrame(cv_interplate).mean(),
                                columns=['CV'])
                remove_met_table = cv_interplate[cv_interplate['CV'] > cutoff]
                metabolite_list.extend(list(remove_met_table.index))
                self._print_metabolites_removed(remove_met_table, i)
                #Remove metabolites in data and pool
                self._remove_metabolites(remove_met_table, i)
            metabolite_list = set(metabolite_list)
            n_metabolites   = len(metabolite_list)
            print('There were ' + 
                  str(n_metabolites) + 
                  ' metabolites removed in the p180 platform across cohorts\n')

    def remove_metabolites_icc(self,
                               cutoff: float = 0.65,
                               by_cohort: bool = True):
        '''
        Compute the intra-class correlation among duplicates or triplicates
        for each metabolite and removes metabolites with ICC lower than cutoff.
        Can only be applied to p180 platform

        Parameters
        ----------
        cutoff: float
            ICC metabolite removal cutoff.
        by_cohort: bool
            Whether to analyze by cohort, or by merging metabolites
            across cohorts (only applicable in p180 platform)

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with metabolites removed due to low ICC.
        pool: pd.Dataframe
            Dataframe with metabolites removed due to low ICC.
        '''
        if self.platform == 'nmr':
            print('Cannot compute the ICC in the NMR platform')
        else:
            print('-----Removing metabolites with ICC values lower than ' +
                  str(cutoff) + '-----')
            metabolite_list = []
            if by_cohort:
                indices = range(len(self.data))
            else:
                indices = [[0,2], [1,3]]
            for i in indices:
                if by_cohort:
                    data = self.data[i]
                else:
                    data = pd.concat([self.data[i[0]],
                                      self.data[i[1]]])
                duplicates_ID  = data.index[\
                                 data.index.duplicated()].unique()
                duplicates_dat = data.loc[duplicates_ID]
                
                raters = []
                for j in duplicates_ID:
                    n_duplicates = len(duplicates_dat.loc[j])
                    for k in range(n_duplicates):
                        raters.append(k+1)
                
                iccs = []
                for met in duplicates_dat.columns:
                    icc_dat = pd.DataFrame()
                    icc_dat['raters']  = raters
                    icc_dat['value']   = list(duplicates_dat[met])
                    icc_dat['targets'] = list(duplicates_dat.index)
                    icc_results = pg.intraclass_corr(icc_dat, 
                                                     targets='targets',
                                                     raters='raters',
                                                     ratings='value',
                                                     nan_policy='omit')
                    iccs.append(icc_results.iloc[2,2])
                
                icc_values = pd.DataFrame(index=duplicates_dat.columns)
                icc_values['ICC'] = iccs
                remove_met_table  = icc_values[icc_values['ICC'] < cutoff]
                metabolite_list.extend(list(remove_met_table.index))
    
                # Print and remove metabolites
                self._print_metabolites_removed(remove_met_table, i)
                self._remove_metabolites(remove_met_table, i)

            metabolite_list = set(metabolite_list)
            n_metabolites   = len(metabolite_list)
            print('There were ' + 
                  str(n_metabolites) + 
                  ' metabolites removed in the p180 platform across cohorts\n')

    def harmonize_metabolites(self):
        '''
        Remove metabolites that have been removed in the other cohort.
        Can only be applied to p180 platform

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with the same metabolites kept across cohorts,
            for the same type
        '''
        if self.platform == 'nmr':
            print('Cannot compute the CV in the NMR platform')
        else:
            print('-----Harmonizing metabolites-----')
            indices = [[0,2],[1,3]]
            for i in indices:
                mets1 = self.data[i[0]].columns
                mets2 = self.data[i[1]].columns
                mets  = [mets1, mets2]
                overlap_mets = mets1.intersection(mets2)
                for l in range(2):
                    remove_mets = mets[l][~mets[l].isin(overlap_mets)]
                    print('We will remove ' +
                          str(len(remove_mets)) + 
                          ' metabolites in ' + 
                          self.cohort[i[l]] + ' ' + 
                          self.type[i[l]]) 
                    self.data[i[l]].drop(remove_mets,
                                         axis=1,
                                         inplace=True)
            print('')

    def impute_metabolites(self):
        '''
        Impute NAs by using 0.5 * LOD score in the p180 platform,
        or by 0 in the nmr platform

        Returns
        ----------
        data: pd.Dataframe
            data with metabolites imputed
        '''
        print('-----Imputing metabolites-----')
        if self.platform == 'p180':
            for i in range(len(self.data)):
                mets_to_impute = self.data[i].\
                                 columns[self.data[i].isna().any()]
                print('We will impute ' +
                      str(len(mets_to_impute)) +
                      ' metabolites in ' +
                      self.cohort[i] + ' ' + 
                      self.type[i])
                for j in mets_to_impute:
                    list_of_rows = self.data[i].loc[\
                                   self.data[i][j].isna()].index
                    for m in list_of_rows:
                        barcode = pd.Series(self.data_meta[i].\
                                            loc[m,'Plate.Bar.Code'])
                        if barcode.size > 1:
                            vals = []
                            for bars in barcode:
                                vals.append(self.lod[i].loc[\
                                     self.lod[i]['Plate.Bar.Code'] == bars, j])
                            self.data[i].at[m, j] = np.mean(vals) * 0.5
                        elif barcode.size == 1:
                            if barcode.isin(self.lod[i]['Plate.Bar.Code']).all():
                                self.data[i].at[m, j] = self.lod[i].loc[\
                                    self.lod[i]['Plate.Bar.Code'].\
                                        isin(barcode), j] * 0.5
                            else:
                                print('The barcode ' +
                                      str(barcode) + 
                                      ' is not in the LOD database for ' +
                                      self.cohort[i] + ' ' + 
                                      self.type[i])
                        else:
                            print('There is something weird here')
        elif self.platform == 'nmr':
            mets_to_impute = self.data.\
                             columns[self.data.isna().any()]
            print('We will impute ' +
                  str(len(mets_to_impute)) + 
                  ' metabolites in the nmr platform')
            self.data[self.data.isna()] = 0
        print('')

    def transform_metabolites_log2(self):
        '''
        Transform metabolite concentration values to log2 values.
        Add a constant of 1 before log transformation.

        Returns
        ----------
        data: pd.Dataframe
            data with metabolite values log2 transformed
        '''
        print('-----Log2 transform-----\n')
        if self.platform == 'p180':
            for i in range(len(self.data)):
                self.data[i] = np.log2(self.data[i] + 1)
        elif self.platform == 'nmr':
            self.data = np.log2(self.data + 1)

    def transform_metabolites_zscore(self):
        '''
        Apply zscore normalization on metabolite values.
        In p180 platform, metabolites are merged across
        cohorts

        Returns
        ----------
        data: pd.DataFrame
            data values normalized
        '''
        print('-----Z-score normalizing data in the ' +
              self.platform +
              ' platform-----\n')
        if self.platform == 'p180':
            indices = [[0,2], [1,3]]
            for i in indices:
                data = pd.concat([self.data[i[0]],
                                  self.data[i[1]]])
                data = zscore_normalize(data)
                for l in i:
                    self.data[l] = data.loc[self.data[l].index]
        elif self.platform == 'nmr':
            self.data = zscore_normalize(self.data)

    def replace_three_std(self):
        '''
        Replace values higher or lower than 3 standard deviations
        with a value of 3 or -3

        Returns
        ----------
        data: pd.Dataframe
            data with values higher or lower than 3 std replaced with
            3 or -3
        '''
        print('-----Replacing extreme values-----')
        for i in range(len(self.data)):
            count = 0
            #three_std   = self.data[i].std() * 3
            three_std   = 3
            bool_higher = self.data[i] > three_std
            bool_lower  = self.data[i] < -three_std
            for row in range(self.data[i].shape[0]):
                for col in range(self.data[i].shape[1]):
                    if bool_higher.iloc[row, col]:
                        count = count + 1
                        self.data[i].iloc[row, col] = three_std
                    elif bool_lower.iloc[row, col]:
                        count = count + 1
                        self.data[i].iloc[row, col] = -three_std
            print('We replaced ' + 
                  str(count) + 
                  ' values in ' +
                  self.cohort[i] + ' ' + 
                  self.type[i])
        print('')
    
    def _print_metabolites_removed(self, 
                                   remove_met_table: pd.DataFrame = None,
                                   index: int or list = None):
        '''
        Print the how many metabolites will be removed, and the value 
        associated with the decision (either missing, CV, or ICC).

        Parameters
        ----------
        remove_met: Dataframe
            Dataframe with metabolites names and either missing, CV or ICC
            values
        index: int or list
            The index value or values to access the correct dataset analyzed
        '''
        if self.platform == 'p180':
            if type(index) == list:
                suffix = self.type[index[0]]
            elif type(index) == int:
                suffix = self.cohort[index] + \
                         ' ' + \
                         self.type[index]
        elif self.platform == 'nmr':
            suffix = self.platform.upper() + \
                     ' platform'

        if len(remove_met_table) == 0:
            print('None of the metabolites were dropped for '+
                  suffix)
        else:
            print('We will remove the following ' +
                  str(len(remove_met_table)) + 
                  ' metabolites for ' + 
                  suffix + 
                  ':')
            print(remove_met_table)

    def _remove_metabolites(self,
                            remove_met_table: pd.DataFrame = None,
                            index: int or list = None):
        '''
        Remove the list of metabolites from remove_met_table indicated 
        in the index
   
        Parameters
        ----------
        remove_met_table: pd.DataFrame
            Dataframe with metabolites names and either missing, CV or ICC
            values
        index: int or list
            The index value or values to access the correct dataset analyzed

        Returns
        ----------
        data: pd.DataFrame
            Dataframe with metabolites removed.
        pool: pd.DataFrame
            Dataframe with metabolites removed.
        '''
        if self.platform == 'p180':
            if type(index) == list:
                for i in index:
                    self.data[i].drop(remove_met_table.index,
                                      axis=1,
                                      inplace=True)
                    self.pool[i].drop(remove_met_table.index,
                                      axis=1,
                                      inplace=True)
            elif type(index) == int:
                self.data[index].drop(remove_met_table.index,
                                      axis=1,
                                      inplace=True)
                self.pool[index].drop(remove_met_table.index,
                                      axis=1,
                                      inplace=True)
        elif self.platform == 'nmr':
            self.data.drop(remove_met_table.index,
                           axis=1,
                           inplace=True)

    def remove_missing_participants(self,
                                    cutoff: float = 0.4):
        '''
        Remove participants with missing data greater than cutoff

        Parameters
        ----------
        cutoff: float
            Missing data removal cutoff. 

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with participants removed
        '''
        print('-----Removing participants with missing data greater than ' +
              str(cutoff) + '-----')
        if self.platform == 'p180':
            indices = [[0,1], [2,3]]
            for i in indices:
                data = pd.concat([self.data[i[0]],
                                  self.data[i[1]]],
                                  axis=1)
                total_part  = data.shape[1]
                remove_part = data.isna().\
                                   sum(axis=1) / total_part > cutoff
                print('We will remove ' +
                      str(sum(remove_part)) + 
                      ' participants for ' + 
                      self.cohort[i[0]] + 
                      ' cohort')
                for l in i:
                    self.data[l] = self.data[l][~remove_part]
        elif self.platform == 'nmr':
            total_part  = self.data.shape[1]
            remove_part = self.data.isna().\
                               sum(axis=1) / total_part > cutoff
            print('We will remove ' +
                  str(sum(remove_part)) + 
                  ' participants for the nmr platform')
            self.data = self.data[~remove_part]
        print('')

    def remove_non_fasters(self):
        '''
        Remove participants that weren't fasting during data collection

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with participants removed 
        '''
        print('-----Removing non-fasting participants-----')
        fasting_participants = self.fasting[self.fasting == 1].index
        if self.platform == 'p180':
            for i in range(len(self.data)):
                keep_participants = self.data[i].index.isin(\
                                    fasting_participants)
                print('We will remove '+
                      str(sum(~keep_participants))+
                      ' participants in '+
                      str(self.cohort[i]) + ' ' +
                      str(self.type[i]))
                self.data[i]      = self.data[i].loc[keep_participants]
        elif self.platform == 'nmr':
            keep_participants = self.data.index.isin(\
                                    fasting_participants)
            print('We will remove '+
                  str(sum(~keep_participants))+
                  ' participants in the nmr platform')
            self.data = self.data.loc[keep_participants]
        print('')

    def remove_multivariate_outliers(self):
        '''
        Remove multivariate outliers by computing the mahalanobis distance
        between participants

        Returns
        ----------
        data: pd.Dataframe
            data with multivariate outliers removed
        '''
        print('-----Removing multivariate outliers-----')
        participant_list = []
        indices = [[0,1],[2,3]]
        for i in indices:
            data = pd.concat([self.data[i[0]],
                              self.data[i[1]]],
                              axis=1)
            cov_mat = data.cov()
            p2      = data.mean()
            cov_mat_pm1 = np.linalg.matrix_power(cov_mat, -1)
            distances = []
            for l, val in enumerate(data.to_numpy()):
                p1 = val
                distance = (p1-p2).T.dot(cov_mat_pm1).dot(p1-p2)
                distances.append(distance)
            distances = np.array(distances)
            cutoff    = stats.chi2.ppf(0.999, data.shape[1]-1)
            n_to_remove = (distances > cutoff ).sum()
            participant_list.extend(list(data.\
                                    index[(distances > cutoff)]))
            print('We will remove ' + 
                  str(n_to_remove) + 
                  ' participants from ' + 
                  self.cohort[i[0]])
            self.data[i[0]]      = self.data[i[0]].\
                                        loc[~(distances > cutoff)]
            self.data[i[1]]      = self.data[i[1]].\
                                        loc[~(distances > cutoff)]
        participant_list = set(participant_list)
        n_participants   = len(participant_list)
        print('There were ' + 
              str(n_participants) + 
              ' that were removed in the p180 platform across cohorts\n')

    def harmonize_participants(self):
        '''
        Remove participants that have been removed in the other cohort

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with the same participants kept across cohorts,
            for the same type
        '''
        print('-----Harmonizing participants-----')
        indices = [[0,1],[2,3]]
        for i in indices:
            parts1 = self.data[i[0]].index
            parts2 = self.data[i[1]].index
            parts  = [parts1, parts2]
            overlap_parts = parts1.intersection(parts2)
            for l in range(2):
                remove_parts = parts[l][~parts[l].isin(overlap_parts)]
                print('We will remove ' +
                      str(len(remove_parts)) + 
                      ' participants in ' + 
                      self.cohort[i[l]] + ' ' + 
                      self.type[i[l]]) 
                self.data[i[l]].drop(remove_parts,
                                     axis=0,
                                     inplace=True)
        print('')

    def consolidate_replicates(self):
        '''
        Consolidate replicates by estimating the average across replicates

        Returns
        ----------
        data: pd.DataFrame
            Dataframe with biological replicates averaged
        '''
        print('-----Consolidating replicates-----')
        if self.platform == 'p180':
            for i in range(len(self.data)):
                duplicates_ID  = self.data[i].index[\
                                 self.data[i].index.duplicated()].unique()
                print('There are ' + 
                      str(len(duplicates_ID))+
                      ' replicated IDs in ' +
                      self.cohort[i] + ' ' + 
                      self.type[i])
                for j in duplicates_ID:
                    consolidated = list(self.data[i].loc[j].mean())
                    self.data[i].drop(j,
                                      axis='index',
                                      inplace=True)
                    self.data[i].loc[j] = consolidated
                self.data[i].sort_index(inplace=True)
        elif self.platform == 'nmr':
            duplicated_IDs = self.data.index[\
                             self.data.index.duplicated()].\
                                  unique()
            print('There are ' + 
                  str(len(duplicated_IDs))+
                  ' replicated IDs in the nmr platform')
            for j in duplicated_IDs:
                consolidated = list(self.data.loc[j].mean())
                self.data.drop(j,
                               axis='index',
                               inplace=True)
                self.data.loc[j] = consolidated
            self.data.sort_index(inplace=True)
        print('')
    
    def compute_cross_plate_correction(self):
        '''
        Computes the cross-plate correction for each metabolite,
        by estimating the average metabolite value across sample pools
        and within plates to generate a per plate correction

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with metabolites corrected 
        '''
        print('-----Computing cross-plate correction-----\n')
        for i in range(len(self.data)):
            global_qc_average = self.pool[i].mean()
            plates_ID = self.pool_meta[i]['Plate.Bar.Code'].unique()
            for j in range(len(plates_ID)):
                # Estimate the correction
                q = self.pool[i][\
                    self.pool_meta[i]['Plate.Bar.Code'] == plates_ID[j]].mean()
                correction = q / global_qc_average

                # Apply the correction
                self.data[i][\
                    self.data_meta[i]['Plate.Bar.Code'] == \
                        plates_ID[j]] = \
                    self.data[i][\
                    self.data_meta[i]['Plate.Bar.Code'] == \
                        plates_ID[j]] / correction
    
    def residuals_from_meds(self,
                            meds):
        '''
        Replace data values with residuals from regression on medication
        dataset.
        Meds dataset is a binary dataframe (presence/absence) of medication
        classes

        Parameters
        ----------
        meds: pd.DataFrame
            medication dataframe
        '''
        print('-----Replacing values with residuals-----\n')
        simplefilter(action='ignore', 
                     category=pd.errors.PerformanceWarning)
        if self.platform == 'p180':
            for i in range(len(self.data)):
                residuals = _get_residuals(self.data[i],
                                           meds.data)
                self.data[i] = residuals
        elif self.platform == 'nmr':
            residuals = _get_residuals(self.data,
                                       meds.data)
            self.data = residuals
    
    def save_files(self):
        '''
        Save cleaned and QCed files to a csv
        '''
        print('-----Saving ' +
              self.platform +
              ' file to .csv-----')
        respath = '../results/'
        if self.platform == 'p180':
            dat1 = pd.concat([self.data[0],
                              self.data[1]],
                              axis=1)
            dat2 = pd.concat([self.data[2],
                              self.data[3]],
                              axis=1)
            final_dat = pd.concat([dat1, dat2])
            name = 'p180_cleaned.csv'
            n_participants = self.data[0].shape[0] + \
                             self.data[3].shape[0]
            n_metabolites  = self.data[0].shape[1] + \
                             self.data[1].shape[1]
        elif self.platform == 'nmr':
            final_dat = self.data
            name = 'nmr_cleaned.csv'
            n_participants = self.data.shape[0]
            n_metabolites  = self.data.shape[1]
        print('Saving ' +
              str(n_participants) +
              ' participants, and ' +
              str(n_metabolites) +
              ' metabolites\n')
        final_dat.to_csv(respath + name)

class QT_pad:
    '''
    QT-PAD dataset class
    '''
    def __init__(self):
        '''
        Initiate the class

        Attributes
        ----------
        data: pd.Dataframe
            qt-pad data
        phenotypes: list
            phenotype column names
        covariates: list
            covariate column names
        extra_vars: list
            extra variables to show in summary
        diagnosis: str
            diagnosis column name
        '''
        print('-----Loading qt-pad data-----')
        self.phenotypes = ['Ventricles',
                           'Hippocampus',
                           'WholeBrain',
                           'Entorhinal',
                           'Fusiform',
                           'MidTemp']
        self.covariates = ['AGE',
                           'PTEDUCAT',
                           'APOE4',
                           'PTGENDER']
        self.extra_vars = ['PTETHCAT',
                           'PTRACCAT']
        self.diagnosis = 'DX.bl'
        qt_path = '../data/ADNI_adnimerge_20170629_QT-freeze.csv'

        dat = pd.read_csv(qt_path).\
                 set_index(['RID','VISCODE'])

        dat = _keep_baseline(dat)

        keep_columns = [self.diagnosis] +\
                        self.covariates +\
                        self.extra_vars +\
                        self.phenotypes +\
                       ['ICV', 'ORIGPROT'] 
                       
        dat = dat.reset_index().set_index('RID')
        dat = dat.loc[:,keep_columns]
        n_participants = dat.shape[0]
        
        print('There are ' +
              str(n_participants) + 
              ' with baseline data')

        self.data = dat.dropna(axis=0)
        n_participants = self.data.shape[0]
        print('There are ' +
              str(n_participants) + 
              ' with complete data\n')
        
        for voi in self.phenotypes:
            voi_corrected  = self.data[voi] / self.data['ICV']
            self.data[voi] = voi_corrected

        # Replace diagnosis categories
        to_replace = {'SMC': 'CN', 
                      'EMCI': 'MCI',
                      'LMCI': 'MCI'}
        self.data[self.diagnosis].replace(to_replace,
                                          inplace=True)
    
    def print_summary(self):
        '''
        Print a summary stats from qtpad data
        '''        
        dat = self.data.copy()
        dat.loc[:,'COUNT'] = 1

        table = pd.pivot_table(dat,
                               index=['ORIGPROT',
                                      'PTGENDER',
                                      'APOE4'],
                               values=['PTEDUCAT',
                                       'AGE',
                                       'COUNT'],
                               aggfunc={'PTEDUCAT': [np.mean, np.std],
                                        'AGE': [np.mean, np.std],
                                        'COUNT': np.sum}).round(1)
        
        print(table)
        print(dat['APOE4'].value_counts())
        print(dat['PTGENDER'].value_counts())
        print(dat['PTRACCAT'].value_counts())
    
    def print_summary_PLS(self):
        '''
        After running a PLS_DA, print a summary with the covariate and 
        phenotype information of the extremes of the PLS scores
        '''
        dat = self.data.copy()
        dat.loc[:,'COUNT'] = 1
        scores = self.scores.copy()

        cols = ['Component 1',
                'Component 2']
        vals = [-1, 1]
        main_print = 'Descriptive statistics for extreme '
        for c in cols:
            for v in vals:
                if c == 'Component 1':
                    threshold = v * 2
                else:
                    threshold = v
                if threshold < 0:
                    extreme_samples = scores[c] < threshold
                    extreme_to_print = 'negative values in '
                elif threshold > 0:
                    extreme_samples = scores[c] > threshold
                    extreme_to_print = 'positive values in '

                table = pd.pivot_table(dat[extreme_samples],
                                       index=['PTGENDER'],
                                       values=['PTEDUCAT',
                                               'AGE',
                                               'Ventricles',
                                               'Hippocampus',
                                               'WholeBrain',
                                               'Entorhinal',
                                               'Fusiform',
                                               'MidTemp',
                                               'COUNT'],
                                       aggfunc={'PTEDUCAT': [np.mean],
                                                'AGE': [np.mean],
                                                'Ventricles': [np.mean],
                                                'Hippocampus': [np.mean],
                                                'WholeBrain': [np.mean],
                                                'Entorhinal': [np.mean],
                                                'Fusiform': [np.mean],
                                                'MidTemp': [np.mean],
                                                'COUNT': np.sum}).round(3)
                
                print(main_print +
                      extreme_to_print +
                      c)
                print(table)
                print('')
        
    def PLS_DA(self,
               n_components:int=2):
        '''
        Run a Partial Least Squares Regression where the
        outcome is the diagnosis, and the predictors are
        the phenotypes

        Parameters
        ----------
        n_components: int
            number of components to estimate

        Returns
        ---------
        scores: pd.DataFrame
            PLS scores
        x_weights: pd.DataFrame
            Weights of X
        y_weights: pd.DataFrame
            Weights of Y
        vips: np.array
            Variable importance in projection
        x_variance_explained: np.array
            Proportion of variance explained by components
        '''
        print('-----Running PLS-DA-----\n')
        dat = self.data[self.phenotypes].\
                   apply(stats.zscore,
                         nan_policy='omit')

        encoder = OneHotEncoder(sparse=False)
        Y = pd.DataFrame(encoder.fit_transform(self.data[ [self.diagnosis] ]))
        Y.columns = encoder.get_feature_names([self.diagnosis])
        
        X = dat.reset_index()[self.phenotypes]
        x_variance = np.sum(np.var(X,
                                   axis = 0))

        PLSDA = PLSRegression(n_components = n_components).\
                              fit(X=X,
                                  Y=Y)
        vips = _calculate_vips(PLSDA)
        self.vips = vips

        x_scores_variance = np.var(PLSDA.x_scores_,
                                   axis = 0) 
                                   
        self.x_variance_explained = x_scores_variance / x_variance
        
        self.scores = pd.DataFrame(PLSDA.x_scores_)
        self.scores.index = self.data.index
        self.scores.rename(columns={0:'Component 1',
                                    1:'Component 2'},
                           inplace=True)
        self.x_weights = pd.DataFrame(PLSDA.x_weights_)
        self.x_weights.index = X.columns
        self.x_weights.rename(columns={0:'Weight 1',
                                       1:'Weight 2'},
                              inplace=True)
        self.y_weights = pd.DataFrame(PLSDA.y_weights_)
        y_weight_columns = ['AD',
                            'CN',
                            'MCI']
        self.y_weights.index = y_weight_columns
        self.y_weights.rename(columns={0:'Weight 1',
                                       1:'Weight 2'},
                              inplace=True)        

    def save_PLS(self,
                 savepath):
        '''
        Save scores from PLS into a csv

        Parameters
        ----------
        savepath: str
            path to save the file
        '''
        print('-----Saving PLS scores to csv-----\n')
        name = 'pheno_pls_components.csv'
        self.scores.to_csv(savepath + name)

class Meds:
    '''
    Medication class
    '''
    def __init__(self):
        '''
        Initiate the class

        Attributes
        ----------
        data: pd.Dataframe
            medication data
        '''
        med_path = '../data/ADMCPATIENTDRUGCLASSES_20170512.csv'

        dat = pd.read_csv(med_path).\
                 set_index(['RID','VISCODE2','Phase'])
        dat.drop('NA',
                 axis=1,
                 inplace=True)

        dat = _keep_baseline(dat)

        self.data = dat
        
    def transform_to_binary(self):
        '''
        Transform the dataset into a binary 0/1 matrix

        Returns
        ----------
        data: pd.Dataframe
            binary matrix
        '''
        # Replacing not NA values
        self.data[self.data.notna()] = 1

        # Replacing NA values
        self.data = self.data.replace(np.nan, 
                                      0)

def zscore_normalize(data,
                     qtpad=None):
    '''
    Apply z-score normalization to data values to mean center 
    and unit variance, by substracting whole column by mean, 
    and dividing by the standard deviation
    Optionally stratify by sex using the qtpad dataset

    Parameters
    ----------
    data: pd.Dataframe
        Dataframe to normalize and RID as index
    qtpad: None or QT_pad
        QT_pad to stratify transformation by sex.
        Participants not in qtpad are removed.

    Returns
    ----------
    normalized_data: pd.Dataframe
        data values normalized
    '''
    normalized_data = []
    if type(data) == list:
        for i in range(len(data)):
            if qtpad is not None:
                dat = pd.merge(qtpad.data['PTGENDER'],
                               data[i],
                               on='RID')
                males_bool   = dat['PTGENDER'] == 'Male'
                females_bool = dat['PTGENDER'] == 'Female'
                males_dat    = dat.drop(['PTGENDER'],
                                        axis=1)[males_bool].\
                                   apply(stats.zscore,
                                         nan_policy='omit')
                females_dat  = dat.drop(['PTGENDER'],
                                        axis=1)[females_bool].\
                                   apply(stats.zscore,
                                         nan_policy='omit')
                final_dat    = pd.concat([females_dat,
                                          males_dat])
            else:
                final_dat = data[i].apply(stats.zscore,
                                          nan_policy='omit')

            normalized_data.append(final_dat)
    else:
        if qtpad is not None:
            dat = pd.merge(qtpad.data['PTGENDER'],
                           data,
                           on='RID')
            males_bool   = dat['PTGENDER'] == 'Male'
            females_bool = dat['PTGENDER'] == 'Female'
            males_dat    = dat.drop(['PTGENDER'],
                                    axis=1)[males_bool].\
                               apply(stats.zscore,
                                     nan_policy='omit')
            females_dat  = dat.drop(['PTGENDER'],
                                    axis=1)[females_bool].\
                               apply(stats.zscore,
                                     nan_policy='omit')
            final_dat    = pd.concat([females_dat,
                                      males_dat])
        else:
            final_dat = data.apply(stats.zscore,
                                   nan_policy='omit')
        normalized_data = final_dat
    
    return(normalized_data)

def _keep_baseline(data):
    '''
    Remove all rows that don't belong to baseline 'bl'

    Parameters
    ----------
    data: pd.Dataframe
        dataframe with VISCODE or VISCODE2 as index

    Returns
    ----------
    baseline_data: pd.Dataframe
        dataframe with only baseline measurements
    '''
    idx = pd.IndexSlice
    baseline_data = data.loc[idx[:, 'bl'], :]

    return(baseline_data)

def _estimate_delete_missingness(data: pd.DataFrame,
                                 cutoff: float = 0.2):
    '''
    Calculate missingness per column and deletes the corresponding
    columns
    Parameters
    ----------
    data: pd.DataFrame
        dataframe with metabolite concentration values
    cutoff: float
        Missing data removal cutoff
    '''
    total_met = len(data)
    remove_met_table = pd.DataFrame(\
                       data.isna().sum(axis=0) / total_met,
                       columns=['missing'])
    remove_met_bool  = data.isna().sum(axis=0) / total_met \
                       > cutoff
    remove_met_table = remove_met_table[remove_met_bool]
    return(remove_met_table)

def _get_residuals(outcomes,
                   predictors):
    '''
    Get residuals for each column in the outcome dataframe.
    Both outcomes and predictors need an 'RID' column

    Parameters
    ----------
    outcomes: pd.DataFrame
        dataframe with the outcome variables
    predictors: pd.DataFrame
        dataframe with the predictor variables
    '''
    new_dat = pd.merge(outcomes,
                       predictors,
                       left_on='RID',
                       right_on='RID')
    Y = new_dat.loc[:,outcomes.columns]
    residuals = Y.copy()
    X = new_dat.loc[:,predictors.columns]
    # Remove meds with only zeros
    keep_meds = X.mean() > 0
    X = X.loc[:,keep_meds]
    for y in Y:
        results = sm.OLS(exog = X,
                         endog = Y[y]).fit()
        med_names = list(X.columns)
        n_significants     = sum(results.pvalues < 0.05)
        n_not_significants = sum(results.pvalues > 0.05)
        while n_not_significants > 0:
            drop_med = results.pvalues[results.pvalues > 0.05].\
                               sort_values(ascending=False).\
                               index[0]
            med_names.remove(drop_med)
            if not med_names:
                print('No significant medications in ' + 
                      y + 
                      '\n')
                break
            else:
                results = sm.OLS(exog = X[med_names],
                                 endog = Y[y]).fit()
                n_significants = sum(results.pvalues < 0.05)
                n_not_significants = sum(results.pvalues > 0.05)
        if n_significants > 0:
            print('There are significant medications in ' +
                  y)
            print(med_names)
            print('')
            residuals[y] = results.resid
    return(residuals)

def _calculate_vips(model):
    '''
    Estimates Variable Importance in Projection (VIP) 
    in Partial Least Squares (PLS)
    
    Parameters
    ----------
    model: PLSRegression
        model generated from the PLSRegression function
    
    Returns
    ----------
    vips: np.array
        variable importance in projection for each variable
    '''
    import numpy as np
    t = model.x_scores_
    w = model.x_weights_
    q = model.y_loadings_
    p, h = w.shape
    vips = np.zeros((p,))
    s = np.diag(np.matmul(np.matmul(np.matmul(t.T,t),q.T), q)).\
           reshape(h, -1)
    total_s = np.sum(s)
    for i in range(p):
        weight = np.array([ (w[i,j] / \
                            np.linalg.norm(w[:,j]))**2 for j in range(h) ])
        vips[i] = np.sqrt(p*(np.matmul(s.T, weight))/total_s)
    return(vips)
