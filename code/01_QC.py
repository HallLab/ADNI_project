import pandas as pd 
import clean

#### Read files
p180  = clean.Metabolites('p180')
nmr   = clean.Metabolites('nmr')
qtpad = clean.QT_pad()
meds  = clean.Meds()

#### Summary information
qtpad.print_summary()

#### Remove metabolites
p180.remove_missing_metabolites()
p180.compute_cross_plate_correction() 
p180.remove_metabolites_cv()
p180.remove_metabolites_icc()
p180.harmonize_metabolites()
nmr.remove_missing_metabolites()

#### Remove participants
p180.remove_missing_participants()
p180.consolidate_replicates()
p180.remove_non_fasters()
p180.harmonize_participants()
nmr.remove_missing_participants()
nmr.consolidate_replicates()
nmr.remove_non_fasters()

#### Data transformation
p180.impute_metabolites()
p180.harmonize_metabolites()
nmr.impute_metabolites()
p180.transform_metabolites_log2()
nmr.transform_metabolites_log2()

p180.transform_metabolites_zscore()
nmr.transform_metabolites_zscore()
p180.replace_three_std()
p180.remove_multivariate_outliers()
p180.harmonize_participants()

#### Residuals from medication
meds.transform_to_binary()
p180.residuals_from_meds(meds)
nmr.residuals_from_meds(meds)

#### Phenotype PLS
qtpad.PLS_DA()
qtpad.remove_multivariate_outliers()
qtpad.save_PLS('../results/')

#### Save files
p180.save_files()
nmr.save_files()
