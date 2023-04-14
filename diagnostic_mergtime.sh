#!/bin/bash

#Este programa procesa los datos de simulacion historica para el periodo 1940/1969 
path='/gws/nopw/j04/ncas_generic/users/jmindlin/output/EP_fluxes_yearly_20230413_144836/work/EP_fluxes/EP_fluxes'
pathout='/home/users/tabu/eunpa_lim_project/EP_FLUX_DATA'
#ACCESS-CM2 ACCESS-ESM1-5 AWI-CM-1-1-MR BCC-ESM1 CanESM5-CanOE CanESM5 CanESM5

diagnostics=(momentum_flux_div momentum_fluxes heat_fluxes heat_flux_div)

for model in AWI-ESM-1-1-LR CESM2-WACCM-FV2 GISS-E2-1-G_r1i1p1f1 IITM-ESM MPI-ESM1-2-HR CanESM5 CNRM-CM6-1 GISS-E2-1-G_r1i1p1f2 IPSL-CM6A-LR MPI-ESM1-2-LR CESM2-FV2 CNRM-CM6-1-HR HadGEM3-GC31-LL KACE-1-0-G MRI-ESM2-0 CESM2-WACCM CNRM-ESM2-1 HadGEM3-GC31-MM MPI-ESM-1-2-HAM NESM3 NorESM2-LM NorESM2-MM SAM0-UNICON TaiESM1 UKESM1-0-LL
do
	for diagnostic in ${diagnostics[@]}
	do	
		list=(`ls ${path}/${model}/${diagnostic}/*.nc`)
	        startfile="${list[0]}#"
        	endfile="${list[-1]}#"
        	startname=${startfile##*/}
        	endname=${endfile##*/}
        	startbase=${startname%_*}
        	endbase=${endname%_*}
        	IFS='_' read -r -a array <<< "$startbase"
        	start=${array[-1]}
        	IFS='-' read -r -a startt <<< "$start"
        	startdate=${startt[0]}
        	IFS='_' read -r -a array <<< "$endbase"
        	end=${array[-1]}
        	IFS='-' read -r -a endt <<< "$end"
        	enddate=${endt[-1]}
        	base=${startname%_*_*}
        	mkdir ${path}/${model}/${diagnostic}/all_times
        	cdo mergetime ${path}/${model}/${diagnostic}/*.nc $path/${model}/${diagnostic}/${base}_$startdate-$enddate"_T42.nc"
	        mv $path/${model}/${diagnostic}/*19*19* $path/${model}/${diagnostic}/*20*20* ${path}/${model}/${diagnostic}/all_times
d	
done; done



