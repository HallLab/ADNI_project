import clean
import clarite
import analyze
import pandas as pd

#### Read results
p180_pheno = analyze.ADNI_sex(metabolite_type='p180')
nmr_pheno = analyze.ADNI_sex(metabolite_type='nmr')

#### Run analysis
p180_pheno.stratified_association()
nmr_pheno.stratified_association()



adni_dat = analyze.ADNI_sex(metabolites, 
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
