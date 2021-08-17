import pandas as pd 
import clean

#### Read files
p180  = clean.Metabolites('p180')
nmr   = clean.Metabolites('nmr')

#### Remove metabolites (P180)
p180.remove_missing_metabolites()
p180.compute_cross_plate_correction() 
p180.remove_metabolites_cv()
p180.remove_metabolites_icc()
p180.harmonize_metabolites()
#### Remove metabolites (NMR)
nmr.remove_missing_metabolites()

#### Remove participants (P180)
p180.remove_missing_participants()
p180.consolidate_replicates()
p180.remove_non_fasters()
p180.harmonize_participants()
#### Remove participants (P180)
nmr.remove_missing_participants()
nmr.consolidate_replicates()
nmr.remove_non_fasters()

#### Data transformation
p180.impute_metabolites()
nmr.impute_metabolites()
p180.transform_metabolites_log2()
nmr.transform_metabolites_log2()

p180.data = clean.zscore_normalize(p180.data)
nmr.data  = clean.zscore_normalize(nmr.data)
p180.replace_three_std()
p180.remove_multivariate_outliers()
p180.harmonize_participants()

#### Medication data
meds = clean.Meds()
meds.transform_to_binary()
p180.residuals_from_meds(meds)
nmr.residuals_from_meds(meds)

#### Save files
p180.save_files()
nmr.save_files()
