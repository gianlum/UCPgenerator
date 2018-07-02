#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

import numpy as np
import gdal
import shapefile

import geo2rot

def lidar(veg_path, thr_val, data1, rlat_v, rlon_v, lat_s_1, lon_s_1, \
        LAD_C, OMEGA, LAI_URB, greening):
    # Inizializations
    # Read the dataset
    print('Reading tree dataset')
    veg = gdal.Open(veg_path)
    # Create lat and lon arrays
    lon_res = veg.GetGeoTransform()[1]
    lon_l = veg.GetGeoTransform()[0]
    lon_r = veg.GetGeoTransform()[0] + veg.RasterXSize * lon_res
    lon = np.arange(lon_l, lon_r, lon_res)
    lat_res = veg.GetGeoTransform()[5]
    lat_l = veg.GetGeoTransform()[3] + lat_res
    lat_r = veg.GetGeoTransform()[3] + veg.RasterYSize * lat_res
    lat = np.arange(lat_l, lat_r, lat_res)
    lon_s_2 = np.zeros((len(lat),len(lon)))
    lat_s_2 = np.zeros((len(lat),len(lon)))
    for j in range(0, len(lat)):
        lon_s_2[j,:] = lon

    for j in range(0, len(lon)):
        lat_s_2[:,j] = lat

    # Read band
    data2 = veg.ReadAsArray()
    data2 = data2.astype(float)
    # Remove anomalous values
    data2[data2>20] = 20 # no trees taller than 20 m
    # Calculations
    print('Calculating LAD') # status info
    OMEGA[0,:,:] = 0.5 # more relistic value for now
    # Write LAD into cosmo output
    for j in range(0, len(lon)):
        for k in range(0, len(lat)):
            h_veg_tmp = data2[k,j]
            if h_veg_tmp > 0.:
                # Identify corresponding grid cell
                [lonV,latV] = np.array(geo2rot.g2r(lon_s_2[k,j],lat_s_2[k,j]))
                lon_idx = np.abs(rlon_v - lonV).argmin()
                lat_idx = np.abs(rlat_v - latV).argmin()
                # Calculate LAD from Lidar canopy height
                area_lidar = 1. # m2 mesured in GIS
                LAD_spec = 1. # m2 m-3 (Klingberg et al., 2017)
                area_grid = 278.93 * 278.93 # m2
                h_ug = 5. # m, height urban grida
                h_cbase = 3. # m, height of the canopy base
                # Check if in-canyon vegetation
                lon_urb_idx = np.abs(lon_s_1[0,:] - lon_s_2[k,j]).argmin()
                lat_urb_idx = np.abs(lat_s_1[:,0] - lat_s_2[k,j]).argmin()
                data1_tmp = data1[lat_urb_idx,lon_urb_idx]
                if data1_tmp >= thr_val:
                    # In-canyon vegetation
                    if h_veg_tmp >= h_cbase :
                        LAD_C[0,:,0,lat_idx,lon_idx] += LAD_spec * area_lidar * \
                            (min(h_veg_tmp,h_ug) - h_cbase) / area_grid / h_ug
                    if h_veg_tmp >= 5 :
                        LAD_C[0,:,1,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,10)-5)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 10 :
                        LAD_C[0,:,2,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,15)-10)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 15 :
                        LAD_C[0,:,3,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,20)-15)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 20 :
                        LAD_C[0,:,4,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,25)-20)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 25 :
                        LAD_C[0,:,5,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,30)-25)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 30 :
                        LAD_C[0,:,6,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,35)-30)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 35 :
                        LAD_C[0,:,7,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,40)-35)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 40 :
                        LAD_C[0,:,8,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,45)-40)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 45 :
                        LAD_C[0,:,9,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,50)-45)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 50 :
                        LAD_C[0,:,10,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,55)-50)/h_ug \
                                                        * area_lidar / area_grid
                    if h_veg_tmp >= 55 :
                        LAD_C[0,:,11,lat_idx,lon_idx] += LAD_spec * (min(h_veg_tmp,60)-55)/h_ug \
                                                        * area_lidar / area_grid

                else:
                    # Out-canyon vegetation
                    LAI_URB[0,lat_idx,lon_idx] += LAD_spec * area_lidar * (h_veg_tmp-h_cbase) \
                                                    / area_grid





    if greening == True: # use greening scenario
        LAD_C = LAD_C * (1 + 2.33 * np.exp(-1 * LAD_C/LAD_C.mean()))

    return LAD_C, OMEGA, LAI_URB


