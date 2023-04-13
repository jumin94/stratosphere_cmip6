import glob
import os
import yaml
import numpy as np 

years = np.arange(1979,2015,1)
keys = ['ua'+str(y) for y in years]+['va'+str(y) for y in years]+['ta'+str(y) for y in years]
var = ['ua']*len(years)+['va']*len(years)+['ta']*len(years)
years2 = np.concatenate([years,years])
years3 = np.concatenate([years2,years])

dic = {}
dic['diagnostics'] = {}
dic['diagnostics']['EP_fluxes'] = {}
dic['diagnostics']['EP_fluxes']['variables'] = {key:{'short_name': v,'mip': 'Amon','start_year':int(str(y)),'end_year':int(str(y)),'preprocessor':'general_preproc'} for key,v,y in zip(keys,var,years3)}


#{'short_name': ua,'mip': Amon,'start_year':y,'end_year':y,'preprocessor':general_preproc},'va'+str(y):{'short_name': va,'mip': Amon,'start_year':y,'end_year':y,'preprocessor':general_preproc},'ta'+str(y):{'short_name': ta,'mip': Amon,'start_year':y,'end_year':y,'preprocessor':general_preproc} for y in years}


with open("/home/users/tabu/recipes/recipe_years.yml", "w") as f:
    yaml.dump(dic, f, default_flow_style=False)               

