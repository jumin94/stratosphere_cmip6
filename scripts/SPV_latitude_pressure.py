#!/usr/bin/env python

#!/usr/bin/env python

import xarray as xr
import numpy as np
import os, fnmatch
import sys
import bottleneck as bn #version 1.0.0
import xarray as xr
import numpy as np
import statsmodels.api as sm
import pandas as pd
import json
import os
from esmvaltool.diag_scripts.shared import run_diagnostic, get_cfg, group_metadata
from scipy import signal
import matplotlib.pyplot as plt

def jet_lat_strength(jet_data,lon1=0,lon2=360):
    jet_data = jet_data.sel(time=is_jja(jet_data['time.month']))
    jet_40_80 = jet_data.mean(dim=time).sel(lat=slice(-80,-40)).sel(lon=slice(lon1,lon2)).mean(dim='lon')**2
    lat = jet_40_80.lat
    jet_lat = (jet_40_80*lat).sum(dim='lat')/(jet_40_80).sum(dim='lat')
    jet_strength = jet_data.groupby('time.year').mean(dim='time').sel(lat=slice(jet_lat-5,jet_lat+5)).sel(lon=slice(lon1,lon2)).mean(dim='lon'))
    return jet_lat,jet_strength

def is_jja(month):
    return (month >= 6) & (month <= 8)

def plot_function(data):
    fig = plt.figure()
    plt.plot(data.time,data.lev)

def main(config):
    """Run the diagnostic."""
    cfg=get_cfg(os.path.join(config["run_dir"],"settings.yml"))
    #print(cfg)
    meta = group_metadata(config["input_data"].values(), "alias")
    meta_year = group_metadata(config["input_data"].values(), "end_year")
    #print(f"\n\n\n{meta}")
    for alias, alias_list in meta.items():
        path_out = config["work_dir"]
        os.chdir(path_out)
        os.getcwd()
        os.makedirs(alias,exist_ok=True)
        os.chdir(path_out+'/'+alias)
        os.getcwd()
        os.makedirs('lat_pressure',exist_ok=True)
        print('alias_list',alias_list)
        u = [xr.open_dataset(m["filename"]) for m in alias_list if m["short_name"] == "ua"]
        #Compute centroid
        centroid_lat, jja_strength = jet_lat_strength(u[0])
        lat_pressure_field = u[0].sel(lat==centroid_lat)
        #Save 
        LP = lat_pressure_field.to_netcdf(path_out+'/'+alias+'/time_pressure/time_pressure_'+alias+'_T42.nc')
        print('Finished with year model '+alias)


		
if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)

