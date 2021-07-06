import clean
import clarite
import analyze
import pandas as pd

#### Read results
p180 = pd.read_csv('~/work/ADNI_project/results/p180_cleaned.csv').\
          set_index('RID')

qtpad = clean.QT_pad()
qtpad.keep_baseline()
qtpad.keep_phenotypes()
qtpad.zscore_normalize_phenotypes()

adni_dat = analyze.ADNI_sex(p180, 
                            qtpad)
adni_dat.run_ewas()
adni_dat.meta_analyze()
adni_dat.sex_diff_test()
adni_dat.categorize_sex_diff()
adni_dat.categorize_arnold20()

savepath = '~/work/ADNI_project/results'
adni_dat.save_to_csv(savepath)
