import xarray as xr
import numpy as np
import os, fnmatch

#t = xr.open_dataset('/home/users/co157815/cmip6_stratosphere/ta_day_ACCESS-ESM1-5_historical_r1i1p1f1_gn_20100101-20141231.nc')
#u = xr.open_dataset('/home/users/co157815/cmip6_stratosphere/ua_day_ACCESS-ESM1-5_historical_r1i1p1f1_gn_20100101-20141231.nc')
#v = xr.open_dataset('/home/users/co157815/cmip6_stratosphere/va_day_ACCESS-ESM1-5_historical_r1i1p1f1_gn_20100101-20141231.nc')

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
	divEPp = dXdp(EPp)
	divEPphi = dXdlat(EPphi)
	divEPflux = divEPp + divEPphi
	return divEPp, divEPphi, divEPflux

def open_data(path,model):
	listOfFiles = os.listdir(path)
	pattern = '*'+model[0]+'*'+model[1]+'*.nc'
	for entry in listOfFiles:
		if fnmatch.fnmatch(entry,pattern):
			data = xr.open_dataset(path+'/'+entry)
			success_pattern = entry
		else:
			print('we do not have '+model[0]+' real '+model[1]+' in '+path)
			data = ' '; pattern = ' '
	return data, success_pattern.split("_",8)

def main():
	#list_models = [['ACCESS-CM2', 'r1i1p1f1'],['ACCESS-ESM1-5', 'r1i1p1f1'],['AWI-CM-1-1-MR', 'r1i1p1f1'],['AWI-ESM-1-1-LR', 'r1i1p1f1'], ['BCC-CSM2-MR', 'r1i1p1f1'],['BCC-ESM1', 'r1i1p1f1'], ['CanESM5-CanOE', 'r1i1p2f1'], ['CanESM5', 'r1i1p1f1'], ['CanESM5', 'r1i1p2f1'], ['CAS-ESM2-0', 'r1i1p1f1'], ['CESM2-FV2', 'r1i1p1f1'], ['CESM2', 'r1i1p1f1'], ['CESM2-WACCM-FV2', 'r1i1p1f1'], ['CESM2-WACCM', 'r1i1p1f1'], ['CIESM', 'r1i1p1f1'], ['CMCC-CM2-HR4', 'r1i1p1f1'], ['CMCC-CM2-SR5', 'r1i1p1f1'], ['CMCC-ESM2', 'r1i1p1f1'], ['CNRM-CM6-1', 'r1i1p1f2'], ['CNRM-CM6-1-HR', 'r1i1p1f2'], ['CNRM-ESM2-1', 'r1i1p1f2'], ['E3SM-1-0', 'r1i1p1f1'], ['E3SM-1-1-ECA', 'r1i1p1f1'], ['E3SM-1-1', 'r1i1p1f1'], ['EC-Earth3-AerChem', 'r1i1p1f1'], ['EC-Earth3-CC', 'r1i1p1f1'], ['EC-Earth3', 'r1i1p1f1'], ['EC-Earth3-Veg', 'r1i1p1f1'], ['EC-Earth3-Veg-LR', 'r1i1p1f1'], ['FGOALS-f3-L', 'r1i1p1f1'], ['FGOALS-g3', 'r1i1p1f1'], ['FIO-ESM-2-0', 'r1i1p1f1'], ['GFDL-ESM4', 'r1i1p1f1'], ['GISS-E2-1-G-CC', 'r1i1p1f1'],['GISS-E2-1-G', 'r1i1p1f1'],['GISS-E2-1-G', 'r1i1p1f2'],['GISS-E2-1-H', 'r1i1p1f1'],['GISS-E2-1-H', 'r1i1p1f2'],['GISS-E2-2-G', 'r1i1p1f1'],['GISS-E2-2-H', 'r1i1p1f1'],['HadGEM3-GC31-LL', 'r1i1p1f3'],['HadGEM3-GC31-MM', 'r1i1p1f3'],['ICON-ESM-LR', 'r1i1p1f1'],['IITM-ESM', 'r1i1p1f1'],['INM-CM5-0', 'r1i1p1f1'],['IPSL-CM5A2-INCA', 'r1i1p1f1'],['IPSL-CM6A-LR', 'r1i1p1f1'],['IPSL-CM6A-LR-INCA', 'r1i1p1f1'],['KACE-1-0-G', 'r1i1p1f1'],['KIOST-ESM', 'r1i1p1f1'],['MIROC6', 'r1i1p1f1'],['MIROC-ES2L', 'r1i1p1f2'],['MPI-ESM-1-2-HAM', 'r1i1p1f1'],['MPI-ESM1-2-HR', 'r1i1p1f1'],['MPI-ESM1-2-LR', 'r1i1p1f1'],['MRI-ESM2-0', 'r1i1p1f1'],['NESM3', 'r1i1p1f1'],['NorCPM1', 'r1i1p1f1'],['NorESM2-LM', 'r1i1p1f1'],['NorESM2-MM', 'r1i1p1f1'],['SAM0-UNICON', 'r1i1p1f1'],['TaiESM1', 'r1i1p1f1'],['UKESM1-0-LL', 'r1i1p1f2'],['UKESM1-1-LL', 'r1i1p1f2']]
	list_models = [['ACCESS-ESM1-5','r1i1p1f1']]
	path = '/home/users/co157815/cmip6_stratosphere'
	path_out = '/home/users/co157815/cmip6_stratosphere'
	for model in list_models:
		#Open data
		t,list_name = open_data(path+'/ta',model)
		v,list_name = open_data(path+'/va',model)
		u,list_name = open_data(path+'/ua',model)
		#Compute
		EPp, EPphi = EP_fluxes(u,v,t)
		div_EPp, div_EPphi, div_EPtot = EP_flux_divergence(EPp,EPphi)
		EPp = EPp.to_dataset(name='EPp')
		EPphi = EPphi.to_dataset(name='EPphi')
		div_EPp = div_EPp.to_dataset(name='div_EPp')
		div_EPphi = div_EPphi.to_dataset(name='div_EPphi')
		div_EPtot = div_EPtot.to_dataset(name='div_EPtot')
		#Save
		EPp = EPp.to_netcdf(path_out+'/heat_fluxes/EP_heat_flux_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6]+'_T42.nc')
		EPphi = EPphi.to_netcdf(path_out+'/momentum_fluxes/EP_momentum_flux_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6]+'_T42.nc')
		div_EPp = div_EPp.to_netcdf(path_out+'/heat_flux_div/EP_heat_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6]+'_T42.nc')
		div_EPphi = div_EPphi.to_netcdf(path_out+'/momentum_flux_div/EP_momentum_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6]+'_T42.nc')
		div_EPtot = div_EPtot.to_netcdf(path_out+'/total_flux_div/EP_total_flux_divergence_'+list_name[1]+'_'+list_name[2]+'_'+list_name[3]+'_'+list_name[4]+'_'+list_name[6]+'_T42.nc')
		print('Finished with model '+model[0])

if __name__ == '__main__':
	main()
