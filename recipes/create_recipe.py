import glob
import os
import yaml
import numpy as np 

years = np.arange(1979,2015,1)

dic = {}
dic['diagnostics'] = {}
dic['diagnostics']['EP_fluxes'] = {}
dic['diagnostics']['EP_fluxes']['variables'] = {'ua'+str(y):{'short_name': ua,'mip': Amon,'start_year':y,'end_year':y,'preprocessor':general_preproc},'va'+str(y):{'short_name': va,'mip': Amon,'start_year':y,'end_year':y,'preprocessor':general_preproc},'ta'+str(y):{'short_name': ta,'mip': Amon,'start_year':y,'end_year':y,'preprocessor':general_preproc} for y in years}


with open("/home/users/tabu/recipes/recipe_years.yml", "w") as f:
    yaml.dump(dic, f, default_flow_style=False)               

