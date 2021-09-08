import clean
import clarite
import analyze
import pandas as pd

#### Loop to read, normalize, analyze and save the results
platforms = ['p180',
             'p180',
             'nmr',
             'nmr']
modules   = [False,
             True,
             False,
             True]

savepath = '../results'
for i in range(0,4):
    if modules:
        end_line = 'with metabolite modules'
    else:
        end_line = 'with single metabolites'
    print('-----Analyzing ' + 
          platforms[i] + 
          ' platform ' +
          end_line)
    dat = analyze.ADNI(metabolite_type=platforms[i],
                       modules=modules[i],
                       phenotypes_pls=True)
    dat.normalize_data()
    dat.stratified_association()
    dat.save_to_csv(savepath)
