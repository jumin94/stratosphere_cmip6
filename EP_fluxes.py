#!/usr/bin/env python

import xarray as xr
import numpy as np
import os, fnmatch
import sys

def tita(T):
        plev = T.plev
        PT = T.ta*(plev/1000)**(-0.286)
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

def dtitadp(tita):
        a = 6371000 #radius in meters
        tita = tita.mean(dim='lon') # less work if first take zonal mean
        dtitadlogp = tita.copy()
        plev = tita.plev.values
        log_p = np.log(plev)
        tita = tita.values.reshape(-1,tita.shape[0],tita.shape[1],tita.shape[2])
        dtitadlogp_arr = c_diff(tita,log_p,2)
        dtitadlogp.values = dtitadlogp_arr[0,:,:,:]
        dtitadp = dtitadlogp/dtitadlogp.plev
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
        asinlat = a*np.sin(lat)
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
        dPTdp= dtitadp(PT)
        lat = U.lat
        a = 6371000 #radius in meters
        coslat = np.cos(lat*3.14/180)
        omega = 7.2921e-5
        f = 2*omega*np.sin(lat)
        U_anom = zonal_anom(U); V_anom = zonal_anom(V); PT_anom = zonal_anom(PT);
        A = zonal_mean(U_anom.ua*V_anom.va)
        B = zonal_mean(PT_anom*V_anom.va)
        EPp = (1/dPTdp)*f*a*B*coslat
        EPphi = -a*A*coslat*coslat
        return EPp, EPphi

def EP_flux_divergence(EPp,EPphi):
	print(EPp)
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

def main(argv):
	path_ua = argv[2]; path_va = argv[3];path_ta = argv[4]
	file_name = argv[1].split("/",14)[-1]
	list_name = file_name.split("_",8)
	model = list_name[2]
	#path_ua = "/".join(path.split("/",13)[:11])+'/ua/'+"/".join(path.split("/",13)[12:])
	path_out = '/gws/nopw/j04/ncas_generic/users/jmindlin/cmip6/historical_EP_fluxes'
	os.chdir(path_out)
	os.getcwd()
	os.makedirs(model,exist_ok=True)
	os.chdir(path_out+'/'+model)
	os.getcwd()
	os.makedirs('heat_fluxes',exist_ok=True)
	os.makedirs('momentum_fluxes',exist_ok=True)
	os.makedirs('heat_flux_div',exist_ok=True)
	os.makedirs('momentum_flux_div',exist_ok=True)
	os.makedirs('total_flux_div',exist_ok=True)
	#Open data
	u = open_data(path_ua)
	v = open_data(path_va)
	t = open_data(path_ta)
	file_name = argv[1].split("/",14)[-1]
	list_name = file_name.split("_",8)
	print(list_name)
	#Compute
	EPp, EPphi = EP_fluxes(u,v,t)
	div_EPp, div_EPphi, div_EPtot = EP_flux_divergence(EPp,EPphi)
	EPp = EPp.to_dataset(name='EPp')
	EPphi = EPphi.to_dataset(name='EPphi')
	div_EPp = div_EPp.to_dataset(name='div_EPp')
	div_EPphi = div_EPphi.to_dataset(name='div_EPphi')
	div_EPtot = div_EPtot.to_dataset(name='div_EPtot')
	print('EP_heat_flux_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
	#Save
	EPp = EPp.to_netcdf(path_out+'/'+list_name[2]+'/heat_fluxes/EP_heat_flux_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
	EPphi = EPphi.to_netcdf(path_out+'/'+list_name[2]+'/momentum_fluxes/EP_momentum_flux_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
	div_EPp = div_EPp.to_netcdf(path_out+'/'+list_name[2]+'/heat_flux_div/EP_heat_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
	div_EPphi = div_EPphi.to_netcdf(path_out+'/'+list_name[2]+'/momentum_flux_div/EP_momentum_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
	div_EPtot = div_EPtot.to_netcdf(path_out+'/'+list_name[2]+'/total_flux_div/EP_total_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6].split(".",2)[0]+'_T42.nc')
	print('Finished with file ',"/".join(list_name))

		
if __name__ == '__main__':
	main(sys.argv)



