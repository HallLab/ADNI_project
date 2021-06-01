import pandas as pd
import pingouin as pg
import numpy as np
from sklearn import preprocessing
from scipy.stats import chi2

class P180:
    '''
    Biocrates p180 metabolomics platform class
    '''
    def __init__(self):
        '''
        Initiate the class by loading all the datasets

        Attributes
        ----------
        cohort: List of str
            List of the ADNI cohort (ADNI1 or ADNI2GO)
        type: List of str
            List of the type of analysis (UPLC or FIA)
        pool: List of Dataframe
            List of data from sample pools
        data: List of Dataframe
            List of data from participants
        data_meta: List of Dataframe
            List of metadata from the participants data
        pool_meta: List of Dataframe
            List of metadata from the pool samples
        '''
        #### SETTING PATHS AND FILES ####
        p180_path = '~/mah546/default/datasets/ADNI/Biospecimen/Biospecimen_Results/ADMC/Biocrates_p180/'
        lod_path  = p180_path + 'LODvalues/'
        fasting_path = '~/mah546/default/datasets/ADNI/Biospecimen/Biospecimen_Results/BIOMARK.csv'

        file_names = ['ADMCDUKEP180UPLC_01_15_16.csv',
                      'ADMCDUKEP180FIA_01_15_16.csv',
                      'ADMCDUKEP180UPLCADNI2GO.csv',
                      'ADMCDUKEP180FIAADNI2GO.csv']
        
        #### METABOLOMICS DATA ####
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
            dat = pd.read_csv(p180_path + i,
                              na_values=na_values).set_index('RID')
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
        
        #### FASTING DATA ####
        fast_dat = pd.read_csv(fasting_path, 
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
                          axis="index",
                          inplace=True)
            fast_dat.loc[i] = val
        self.fasting = fast_dat.sort_index()

        #### LOD VALUES ####
        for i in range(len(self.data)):
            filename = 'P180' + \
                        str(self.type[i]) + \
                        'LODvalues_' + \
                        str(self.cohort[i]) + \
                        '.csv'
            dat = pd.read_csv(lod_path + filename)
            # In lod value ADNI2GO FIA, the bar code plate
            # needs fixing
            if i == 3:
                barcode = dat['Plate Bar Code']
                barcode = barcode.str.split(' ',
                                            expand = True)[2].\
                                  str.replace(pat='/',
                                              repl='-')
                dat['Plate Bar Code'] = barcode

            # Metabolite names in lod don't match those in
            # the data, replace '-', ':', '(', ')' and ' ' with '.'
            old_columns = dat.columns
            new_columns = old_columns.str.replace(pat='-|:|\(|\)| ', repl='.')
            dat.columns = new_columns
            # Change metabolite name from Met.So to Met.SO
            if self.type[i] == 'UPLC':
                dat.rename(columns={'Met.SO':'Met.So'}, inplace=True)
            self.lod.append(dat)

    def remove_missing_metabolites(self,
                                   cutoff: float = 0.2):
        '''
        Remove metabolites due to missing data greater than cutoff
        
        Parameters
        ----------
        cutoff: float
            Missing data removal cutoff. 

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with metabolites removed due to missingness
        pool: pd.Dataframe
            Dataframe with metabolites removed due to missingness
        '''
        print('-----Removing metabolites with missing data greater than ' +
              str(cutoff) + '-----')
        for i in range(len(self.data)):
            total_met = len(self.data[i])
            remove_met_table = pd.DataFrame(\
                               self.data[i].isna().sum(axis=0) / total_met,
                               columns=['missing'])
            remove_met_bool  = self.data[i].isna().sum(axis=0) / total_met \
                               > cutoff
            remove_met_table = remove_met_table[remove_met_bool]
            self._print_metabolites_removed(remove_met_table, i)
            #Remove metabolites from data and pool
            self._remove_metabolites(remove_met_table, i)

    def remove_metabolites_cv(self,
                              cutoff: float = 0.2):
        '''
        Compute the coefficient of variation among duplicates or triplicates
        for each metabolite and remove metabolites with CV higher than cutoff

        Parameters
        ----------
        cutoff: float
            CV metabolite removal cutoff. 

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with metabolites removed due to high CV.
        pool: pd.Dataframe
            Dataframe with metabolites removed due to high CV.
        '''
        print('-----Removing metabolites with CV values greater than ' +
              str(cutoff) + '-----')
        for i in range(len(self.data)):
            cv_interplate = []
            duplicates_ID = self.data[i].index[\
                            self.data[i].index.duplicated()].unique()
            for j in range(len(duplicates_ID)):
                duplicates_dat = self.data[i].loc[duplicates_ID[j],:]
                cv = duplicates_dat.std() / duplicates_dat.mean()
                cv_interplate.append(cv)
            cv_interplate = pd.DataFrame(pd.DataFrame(cv_interplate).mean(),
                            columns=['CV'])
            remove_met_table = cv_interplate[cv_interplate['CV'] > cutoff]
            self._print_metabolites_removed(remove_met_table, i)
            #Remove metabolites in data and pool
            self._remove_metabolites(remove_met_table, i)

    def remove_metabolites_icc(self,
                               cutoff: float = 0.65):
        '''
        Compute the intra-class correlation among duplicates or triplicates
        for each metabolite and removes metabolites with ICC lower than cutoff

        Parameters
        ----------
        cutoff: float
            ICC metabolite removal cutoff. 

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with metabolites removed due to low ICC.
        pool: pd.Dataframe
            Dataframe with metabolites removed due to low ICC.
        '''
        print('-----Removing metabolites with ICC values lower than ' +
              str(cutoff) + '-----')
        for i in range(len(self.data)):
            duplicates_ID  = self.data[i].index[\
                             self.data[i].index.duplicated()].unique()
            duplicates_dat = self.data[i].loc[duplicates_ID]
            
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

            # Print and remove metabolites
            self._print_metabolites_removed(remove_met_table, i)
            self._remove_metabolites(remove_met_table, i)

    def harmonize_metabolites(self):
        '''
        Remove metabolites that have been removed in the other cohort

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with the same metabolites kept across cohorts,
            for the same type
        '''
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
                      self.type[i[l]] + '\n' ) 
                self.data[i[l]].drop(remove_mets,
                                     axis=1,
                                     inplace=True)

    def impute_metabolites(self):
        '''
        Impute NAs by using 0.5 * LOD score

        Returns
        ----------
        data: pd.Dataframe
            data with metabolites imputed
        '''
        print('-----Imputing metabolites-----')
        for i in range(len(self.data)):
            mets_to_impute = self.data[i].\
                             columns[self.data[i].isna().any()]
            print('We will impute ' +
                  str(len(mets_to_impute)) +
                  ' metabolites in ' +
                  self.cohort[i] + ' ' + 
                  self.type[i] + '\n')
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

    def transform_metabolites_log2(self):
        '''
        Transform metabolites concentration values to log2 values

        Returns
        ----------
        data: pd.Dataframe
            data with metabolite values log2 transformed
        '''
        print('-----Log2 transform-----\n')
        for i in range(len(self.data)):
            self.data[i] = np.log2(self.data[i])

    def scale_metabolites(self):
        '''
        Scale metabolites values to mean center and unit variance,
        by substracting whole column by mean, and dividing by the 
        standard deviation

        Returns
        ----------
        data: pd.Dataframe
            data with scaled values
        '''
        print('-----Scaling metabolites-----\n')
        for i in range(len(self.data)):
            preprocessing.scale(self.data[i],
                                axis=0,
                                with_mean=True,
                                with_std=True,
                                copy=False)
    
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
                  self.type[i] + '\n')

    def _print_metabolites_removed(self, 
                                   remove_met_table: pd.DataFrame = None,
                                   index: int = None):
        '''
        Print the how many metabolites will be removed, and the value 
        associated with the decision (either missing, CV, or ICC).

        Parameters
        ----------
        remove_met: Dataframe
            Dataframe with metabolites names and either missing, CV or ICC
            values
        index: int
            The index value to access the correct dataset analyzed
        '''
        if len(remove_met_table) == 0:
            print('None of the metabolites were dropped for '+
                  self.cohort[index] + ' ' +
                  self.type[index] + '\n')
        else:
            print('We will remove the following ' +
                  str(len(remove_met_table)) + 
                  ' metabolites for ' + 
                  self.cohort[index] + ' ' +
                  self.type[index] + ':')
            print(remove_met_table)
            print('')

    def _remove_metabolites(self,
                            remove_met_table: pd.DataFrame = None,
                            index: int = None):
        '''
        Remove the list of metabolites from remove_met_table indicated 
        in the index
   
        Parameters
        ----------
        remove_met: Dataframe
            Dataframe with metabolites names and either missing, CV or ICC
            values
        index: int
            The index value to access the correct dataset analyzed

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with metabolites removed.
        pool: pd.Dataframe
            Dataframe with metabolites removed.
        '''
        self.data[index].drop(remove_met_table.index, axis=1, inplace=True)
        self.pool[index].drop(remove_met_table.index, axis=1, inplace=True)

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
        for i in range(len(self.data)):
            total_part  = self.data[i].shape[1]
            remove_part = self.data[i].isna().sum(axis=1) / total_part > cutoff
            print('We will remove ' +
                  str(sum(remove_part)) + 
                  ' participants for ' + 
                  self.cohort[i] + ' ' +
                  self.type[i] + '\n')
            self.data[i]      = self.data[i][~remove_part]

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
        for i in range(len(self.data)):
            keep_participants = self.data[i].index.isin(\
                                fasting_participants)
            print('We will remove '+
                  str(sum(~keep_participants))+
                  ' participants in '+
                  str(self.cohort[i]) + ' ' +
                  str(self.type[i]) + '\n')
            self.data[i]      = self.data[i].loc[keep_participants]

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
        for i in range(len(self.data)):
            cov_mat = self.data[i].cov()
            p2      = self.data[i].mean()
            cov_mat_pm1 = np.linalg.matrix_power(cov_mat, -1)
            distances = []
            for l, val in enumerate(self.data[i].to_numpy()):
                p1 = val
                distance = (p1-p2).T.dot(cov_mat_pm1).dot(p1-p2)
                distances.append(distance)
            distances = np.array(distances)
            cutoff    = chi2.ppf(0.99, self.data[i].shape[1])
            n_to_remove = (distances > cutoff ).sum()
            print('We will remove ' + 
                  str(n_to_remove) + 
                  ' participants from ' + 
                  self.cohort[i] + ' ' + 
                  self.type[i] + '\n')
            self.data[i]      = self.data[i].loc[~(distances > cutoff)]

    def harmonize_participants(self):
        '''
        Remove participants that have been removed in the other cohort

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with the same participants kept across cohorts,
            for the same type
        '''
        print('-----Harmonizing participants-----\n')
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
                      self.type[i[l]] + '\n' ) 
                self.data[i[l]].drop(remove_parts,
                                     axis=0,
                                     inplace=True)

    def consolidate_replicates(self):
        '''
        Consolidate replicates by estimating the average across replicates

        Returns
        ----------
        data: pd.Dataframe
            Dataframe with biological replicates averaged
        '''
        print('-----Consolidating replicates-----')
        for i in range(len(self.data)):
            duplicates_ID  = self.data[i].index[\
                             self.data[i].index.duplicated()].unique()
            print('There are ' + 
                  str(len(duplicates_ID))+
                  ' replicated IDs in ' +
                  self.cohort[i] + ' ' + 
                  self.type[i] + '\n')
            for j in duplicates_ID:
                consolidated = list(self.data[i].loc[j].mean())
                self.data[i].drop(j,
                                  axis="index",
                                  inplace=True)
                self.data[i].loc[j] = consolidated
            self.data[i].sort_index(inplace=True)
    
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
