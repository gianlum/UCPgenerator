#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

import numpy as np
import gdal
import copy

from . import geo2rot

def wgs(ufrac_path, threshold, thr_val, rlat_v, rlon_v, FR_URBAN):
    # Inizializations
    FR_URBAN_count = copy.deepcopy(FR_URBAN)
    # Read the dataset
    ufrac = gdal.Open(ufrac_path)
    # Create lat and lon arrays
    lon_res = ufrac.GetGeoTransform()[1]
    lon_l = ufrac.GetGeoTransform()[0]
    lon_r = ufrac.GetGeoTransform()[0] + ufrac.RasterXSize * lon_res
    lon = np.arange(lon_l, lon_r, lon_res)
    lat_res = ufrac.GetGeoTransform()[5]
    lat_l = ufrac.GetGeoTransform()[3] + lat_res
    lat_r = ufrac.GetGeoTransform()[3] + ufrac.RasterYSize * lat_res
    lat = np.arange(lat_l, lat_r, lat_res)
    lon_s_1 = np.zeros((len(lat),len(lon)))
    lat_s_1 = np.zeros((len(lat),len(lon)))
    for j in range(0, len(lat)):
        lon_s_1[j,:] = lon

    for j in range(0, len(lon)):
        lat_s_1[:,j] = lat

    # Read band
    data1 = ufrac.ReadAsArray()
    data1 = data1.astype(float)
    # Write FR_URB into cosmo output
    for j in range(0, len(lon)):
        for k in range(0, len(lat)):
            data1_tmp = data1[k,j]
            # Identify corresponding grid cell
            [lonC,latC] = np.array(geo2rot.g2r(lon_s_1[k,j],lat_s_1[k,j]))
            lon_idx = np.abs(rlon_v - lonC).argmin()
            lat_idx = np.abs(rlat_v - latC).argmin()
            # Calculate urban fraction contribution
            data1_tmp = data1_tmp / 100.
            FR_URBAN_count[0,lat_idx,lon_idx] += 1.
            if threshold == 1:
                if data1_tmp >= thr_val:
                    FR_URBAN[0,lat_idx,lon_idx] += 1.
    
            elif threshold == 0:
                FR_URBAN[0,lat_idx,lon_idx] += data1_tmp



    # Disable division by 0 error
    np.seterr(divide='ignore')
    FR_URBAN = FR_URBAN/FR_URBAN_count

    return FR_URBAN, data1, lat_s_1, lon_s_1


