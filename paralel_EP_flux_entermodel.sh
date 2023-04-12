#!/bin/bash

#Este programa procesa los datos de simulacion historica para el periodo 1940/1969 
path='/badc/cmip6/data/CMIP6/CMIP'
path_out='/gws/nopw/j04/ncas_generic/users/jmindlin/cmip6/historical_EP_fluxes'

modelin=$1

echo "$modelin"
#tmp=${a#*_}   # remove prefix ending in "_"
#b=${tmp%_*}   # remove suffix starting with "_"
#AWI-ESM-1-1-LR BCC-CSM2-MR TaiESM1 FGOALS-f3-L FGOALS-g3 ACCESS-CM2 ACCESS-ESM1-5 CanESM5 CanESM5 CESM2-FV2 CESM2 CNRM-CM6-1 CMCC-ESM2 CMCC-CM2-SR5 CMCC-CM2-HR4 CIESM CESM2-WACCM CESM2-WACCM-FV2 CNRM-CM6-1-HR CNRM-ESM2-1 EC-Earth3-AerChem EC-Earth3-CC EC-Earth3 EC-Earth3-Veg-LR FIO-ESM-2-0 GFDL-ESM4 GISS-E2-1-G-CC GISS-E2-1-G GISS-E2-1-G GISS-E2-1-H GISS-E2-1-H GISS-E2-2-G GISS-E2-2-H HadGEM3-GC31-LL HadGEM3-GC31-MM ICON-ESM-LR IITM-ESM INM-CM5-0 IPSL-CM5A2-INCA IPSL-CM6A-LR IPSL-CM6A-LR-INCA KACE-1-0-G KIOST-ESM MIROC6 MIROC-ES2L MPI-ESM-1-2-HAM MPI-ESM1-2-HR MPI-ESM1-2-LR MRI-ESM2-0 NESM3 NorCPM1 NorESM2-LM NorESM2-MM SAM0-UNICON UKESM1-0-LL UKESM1-1-LL EC-Earth3-Veg
#AWI-CM-1-1-MR CanESM5-CanOE CAS-ESM2-0 CIESM E3SM-1-0 E3SM-1-1-ECA E3SM-1-1

dataList_u=()
dataList_v=()
dataList_t=()

mkdir $path_out/$modelin

for model in $modelin
do
for dir in $path/*/
do   dir=${dir%*/}                                              #GRUPO
    for dir2 in ${dir}/*${model}*/
    do   dir2=${dir2%*/}                                        #MODELO
        for dir3 in ${dir2}/historical/r1i*/
        do   dir3=${dir3%*/}                                    #REALIZACION
            for dir4 in ${dir3}/day/ua/*/
            do   dir4=${dir4%*/}
                for filename in ${dir4}/latest/*.nc             #ARCHIVO
                do
                    name=${filename##*/}
                    base=${name%_*}
                    file=$path_out/"${base}.nc"
		    path_out_out=$path_out/"${model}"
                    if [ -f $filename ]; then
                        echo "$filename exists."
			dataList_u+=($filename)
			echo "${base_temp}.nc"
                    else 
                        echo "$filename does not exist."
                    fi
                done
done; done; done; done; done

for model in $modelin
do
for dir in $path/*/
do   dir=${dir%*/}                                              #GRUPO
    for dir2 in ${dir}/*${model}*/
    do   dir2=${dir2%*/}                                        #MODELO
        for dir3 in ${dir2}/historical/r1i*/
        do   dir3=${dir3%*/}                                    #REALIZACION
            for dir4 in ${dir3}/day/va/*/
            do   dir4=${dir4%*/}
                for filename in ${dir4}/latest/*.nc             #ARCHIVO
                do
                    name=${filename##*/}
                    base=${name%_*}
                    file=$path_out/"${base}.nc"
                    path_out_out=$path_out/"${model}"
                    if [ -f $filename ]; then
                        echo "$filename exists."
                        dataList_v+=($filename)
                        echo "${base_temp}.nc"
                    else
                        echo "$filename does not exist."
                    fi
                done
done; done; done; done; done

for model in $modelin
do
for dir in $path/*/
do   dir=${dir%*/}                                              #GRUPO
    for dir2 in ${dir}/*${model}*/
    do   dir2=${dir2%*/}                                        #MODELO
        for dir3 in ${dir2}/historical/r1i*/
        do   dir3=${dir3%*/}                                    #REALIZACION
            for dir4 in ${dir3}/day/ta/*/
            do   dir4=${dir4%*/}
                for filename in ${dir4}/latest/*.nc             #ARCHIVO
                do
                    name=${filename##*/}
                    base=${name%_*}
                    file=$path_out/"${base}.nc"
                    path_out_out=$path_out/"${model}"
                    if [ -f $filename ]; then
                        echo "$filename exists."
                        dataList_t+=($filename)
                        echo "${base_temp}.nc"
                    else
                        echo "$filename does not exist."
                    fi
                done
done; done; done; done; done

for i in "${dataList_u[@]}"
do
     #cdo remapbil,n32 value $path_out/"temp_T42.nc"
     echo i
done

for i in "${!dataList_u[@]}"
do
     cdo remapbil,n32 "${dataList_u[i]}" $path_out/$modelin/"temp_u_T42.nc"
     cdo remapbil,n32 "${dataList_v[i]}" $path_out/$modelin/"temp_v_T42.nc"
     cdo remapbil,n32 "${dataList_t[i]}" $path_out/$modelin/"temp_t_T42.nc"
     python /home/users/tabu/eunpa_lim_project/EP_fluxes.py "${dataList_u[i]}" $path_out/$modelin/"temp_u_T42.nc" $path_out/$modelin/"temp_v_T42.nc" $path_out/$modelin/"temp_t_T42.nc"
     rm $path_out/$modelin/"temp_u_T42.nc" $path_out/$modelin/"temp_v_T42.nc" $path_out/$modelin/"temp_t_T42.nc"
done

