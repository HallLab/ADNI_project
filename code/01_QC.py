import pandas as pd 
import clean

#### Read files
p180 = clean.P180()

#### Remove metabolites
p180.remove_missing_metabolites()
p180.compute_cross_plate_correction() 
p180.remove_metabolites_cv()
p180.remove_metabolites_icc()
p180.harmonize_metabolites()

#### Remove participants
p180.remove_missing_participants()
p180.consolidate_replicates()
p180.remove_non_fasters()
p180.harmonize_participants()

#### Data transformation
p180.impute_metabolites()
p180.transform_metabolites_log2()
p180.zscore_normalize_metabolites()
p180.replace_three_std()
p180.remove_multivariate_outliers()
p180.harmonize_participants()

#### Save files
p180.save_files()
