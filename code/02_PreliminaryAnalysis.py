import clean
import clarite
import analyze
import pandas as pd

#### Read results
p180 = pd.read_csv('../results/p180_cleaned.csv').\
          set_index('RID')
eigenmet = pd.read_csv('../results/eigenmetabolites.csv')
eigenmet = eigenmet.set_index(p180.index)
p180 = pd.merge(p180,
                eigenmet,
                on='RID')
qtpad = clean.QT_pad()
qtpad.keep_baseline()
qtpad.keep_phenotypes()
qtpad.zscore_normalize_phenotypes()

adni_dat = analyze.ADNI_sex(p180, 
                            qtpad)
adni_dat.zscore_normalize_eigenmetabolites()

# Run analysis on eigenmetabolites
adni_dat.run_ewas(on_modules=True)
adni_dat.meta_analyze()
adni_dat.sex_diff_test()
adni_dat.categorize_sex_diff(on_modules=True)
adni_dat.categorize_arnold20(on_modules=True)

savepath = '../results'
adni_dat.save_to_csv(savepath,
                     modules=True)

# Run analysis on individual metabolites
adni_dat.run_ewas()
adni_dat.meta_analyze()
adni_dat.sex_diff_test()
adni_dat.categorize_sex_diff()
adni_dat.categorize_arnold20()

adni_dat.save_to_csv(savepath)
