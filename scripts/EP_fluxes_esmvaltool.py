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

def tita(T):
        plev = T.plev
        PT = T.ta*(plev/100000)**(-0.286)
        return PT

def c_diff(arr, h, dim, cyclic = False):
        ndim = arr.ndim
        lst = [i for i in range(ndim)]
        lst[dim], lst[0] = lst[0], lst[dim]
        rank = lst
        arr = np.transpose(arr, tuple(rank))
        if ndim == 3:
                shp = (arr.shape[0]-2,1,1)
        elif ndim == 4:
                shp = (arr.shape[0]-2,1,1,1)
        elif ndim == 5:
                shp = (arr.shape[0]-2,1,1,1,1)
        d_arr = np.copy(arr)
        if not cyclic:
                d_arr[0,...] = (arr[1,...]-arr[0,...])/(h[1]-h[0])
                d_arr[-1,...] = (arr[-1,...]-arr[-2,...])/(h[-1]-h[-2])
                d_arr[1:-1,...] = (arr[2:,...]-arr[0:-2,...])/np.reshape(h[2:]-h[0:-2], shp)
        elif cyclic:
                d_arr[0,...] = (arr[1,...]-arr[0,...])/(h[1]-h[-1])
                d_arr[-1,...] = (arr[-1,...]-arr[-2,...])/(h[0]-h[-2])
                d_arr[1:-1,...] = (arr[2:,...]-arr[0:-2,...])/np.reshape(h[2:]-h[0:-2], shp)
        d_arr = np.transpose(d_arr, tuple(rank))
        return d_arr

def c_diff_one(arr, h):
#	print arr.shape
	d_arr = np.copy(arr)
	d_arr[0,...] = (arr[1,...]-arr[0,...])/(h[1]-h[0])
	d_arr[-1,...] = (arr[-1,...]-arr[-2,...])/(h[-1]-h[-2])
	d_arr[1:-1,...] = (arr[2:,...]-arr[0:-2,...])/(np.reshape(h[2:]-h[0:-2],(arr.shape[0]-2,1,1)))
	return d_arr


def dtitadp(tita):
        a = 6371000 #radius in meters
        tita = tita.mean(dim='lon') # less work if first take zonal mean
        dtitadlogp = tita.copy()
        plev = tita.plev.values
        log_p = np.log(plev)
        tita = tita.values.reshape(-1,tita.shape[0],tita.shape[1],tita.shape[2])
        dtitadlogp_arr = c_diff(tita,log_p,2)
        dtitadlogp.values = dtitadlogp_arr[0,:,:,:]
        dtitadp = dtitadlogp/np.log(dtitadlogp.plev)
        return dtitadp

def dXdp(X):
        dxdlogp = X.copy()
        plev = X.plev.values
        log_p = np.log(plev)
        X = X.values.reshape(-1,X.shape[0],X.shape[1],X.shape[2])
        dxdlogp_arr = c_diff(X,log_p,2)
        dxdlogp.values = dxdlogp_arr[0,:,:,:]
        dxdp = dxdlogp/dxdlogp.plev
        return dxdp

def dXdlat(X):
        dxdlat = X.copy()
        lat = X.lat.values
        a = 6371000 #radius in meters
        PI = np.pi
        phi = lat*PI/180.0
        asinlat = a*np.sin(phi)
        X = X.values.reshape(-1,X.shape[0],X.shape[1],X.shape[2])
        dxdasinlat_arr = c_diff(X,asinlat,3)
        dxdlat.values = dxdasinlat_arr[0,:,:,:]
        return dxdlat

def zonal_anom(X):
        return X - X.mean(dim='lon')

def zonal_mean(X):
        return X.mean(dim='lon')

def EP_fluxes(U,V,T):
        PT = tita(T)
        #dPTdp= dtitadp(PT)
        dPTdp = PT.mean(dim='lon').copy()
        theta = T.ta.values*(np.reshape(T.plev.values,(1,len(T.plev.values),1,1))/100000.)**(-0.286)
        theta_zm = bn.nanmean(theta, axis = 3)
        loglevel = np.log(T.plev.values)

        dPTdp_arr = np.transpose(c_diff_one(np.transpose(theta_zm,[1,0,2]), loglevel),[1,0,2])
        dPTdp_arr /= 10000.*np.reshape(T.plev.values,(1,len(T.plev.values),1))
        dPTdp.values = dPTdp_arr

        lon = T.lon.values #dataset.variables['lon'][:]
        lat = U.lat
        
        #constants
        a = 6.37122e06 
        PI = np.pi
        phi = lat*PI/180.0     
        acphi=a*np.cos(phi)       
        asphi=a*np.sin(phi)       
        omega = 7.2921e-5      
        f = 2*omega*np.sin(phi) 
        latfac=acphi*np.cos(phi)
        
        U_anom = zonal_anom(U); V_anom = zonal_anom(V); PT_anom = zonal_anom(PT);
        A = zonal_mean(U_anom.ua*V_anom.va)
        B = zonal_mean(PT_anom*V_anom.va)
        EPp = (1/dPTdp)*f*B*acphi
        EPphi = -A*latfac
        return EPp, EPphi

def EP_flux_divergence(EPp,EPphi):
	#print(EPp)
	divEPp = dXdp(EPp)
	divEPphi = dXdlat(EPphi)
	divEPflux = divEPp + divEPphi
	return divEPp, divEPphi, divEPflux

def get_data_files(path,model):
	listOfFiles = os.listdir(path)
	pattern = '*'+model+'*.nc'
	list_files = []        
	for entry in listOfFiles:
		if fnmatch.fnmatch(entry,pattern):
			data = xr.open_dataset(path+'/'+entry)
			list_files.append(entry)
		else:
			print('we do not have '+model[0]+' real '+model[1]+' in '+path)
			data = ' '; pattern = ' '
	return list_files

def open_data(path):
	return xr.open_dataset(path)

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
        os.makedirs('heat_fluxes',exist_ok=True)
        os.makedirs('momentum_fluxes',exist_ok=True)
        os.makedirs('heat_flux_div',exist_ok=True)
        os.makedirs('momentum_flux_div',exist_ok=True)
        os.makedirs('total_flux_div',exist_ok=True)
        print('alias_list',alias_list)
        #Open data
        for year,year_list in meta_year.items():
            print(year)
            u = [xr.open_dataset(m["filename"]) for m in alias_list if m["short_name"] == "ua" and m["end_year"] == year]
            v = [xr.open_dataset(m["filename"]) for m in alias_list if m["short_name"] == "va" and m["end_year"] == year]
            t = [xr.open_dataset(m["filename"]) for m in alias_list if m["short_name"] == "ta" and m["end_year"] == year]
            file_name = [m["filename"].split("/",12)[-1] for m in alias_list if m["short_name"] == "ta" and m["end_year"] == year]
            list_name = file_name[0].split("_",6)
            print(list_name)
            #Compute
            EPp, EPphi = EP_fluxes(u[0],v[0],t[0])
            div_EPp, div_EPphi, div_EPtot = EP_flux_divergence(EPp,EPphi)
            EPp = EPp.to_dataset(name='EPp')
            EPphi = EPphi.to_dataset(name='EPphi')
            div_EPp = div_EPp.to_dataset(name='div_EPp')
            div_EPphi = div_EPphi.to_dataset(name='div_EPphi')
            div_EPtot = div_EPtot.to_dataset(name='div_EPtot')
            print('EP_heat_flux_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
            #Save
            EPp = EPp.to_netcdf(path_out+'/'+alias+'/heat_fluxes/EP_heat_flux_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
            EPphi = EPphi.to_netcdf(path_out+'/'+alias+'/momentum_fluxes/EP_momentum_flux_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
            div_EPp = div_EPp.to_netcdf(path_out+'/'+alias+'/heat_flux_div/EP_heat_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
            div_EPphi = div_EPphi.to_netcdf(path_out+'/'+alias+'/momentum_flux_div/EP_momentum_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
            div_EPtot = div_EPtot.to_netcdf(path_out+'/'+alias+'/total_flux_div/EP_total_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
            print('Finished with year '+str(year)+'model '+alias)


		
if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
