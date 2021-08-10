import pandas as pd 
import clean

#### Read files
p180  = clean.P180()
qtpad = clean.QT_pad()
nmr   = clean.NMR()

#### Remove metabolites (P180)
p180.remove_missing_metabolites()
p180.compute_cross_plate_correction() 
p180.remove_metabolites_cv()
p180.remove_metabolites_icc()
p180.harmonize_metabolites()

#### Remove participants (P180)
p180.remove_missing_participants()
p180.consolidate_replicates()
p180.remove_non_fasters()
p180.harmonize_participants()

#### Data transformation
p180.impute_metabolites()
p180.transform_metabolites_log2(qtpad)
p180.data = clean.zscore_normalize(p180.data, 
                                   qtpad)
nmr.average_replicates()
nmr.metabolites = clean.zscore_normalize(nmr.metabolites, 
                                         qtpad)
p180.replace_three_std()
p180.remove_multivariate_outliers()
p180.harmonize_participants()

#### Medication data
meds = clean.Meds()
meds.transform_to_binary()
p180.residuals_from_meds(meds)

#### Save files
p180.save_files()
nmr.save_files()
