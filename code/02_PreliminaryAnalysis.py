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
    dat = analyze.ADNI(metabolite_type=platforms[i],
                       modules=modules[i])
    dat.normalize_data()
    dat.stratified_association()
    dat.save_to_csv(savepath)
