#!python

import os
import glob
import gdal
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
import numpy.ma as ma

name = 'AHF_Basel_domain'
jday = 173 # Julian Day
ndays = 4
sp = 48 # average vehicle speed
frv = 0.2 # fraction of vehicles on the rood
ext = '.asc.tif'

data_path = '/project/mugi/nas/PAPER2/CCLM-DCEP-Tree/ahf/tif/'
cosmo_path = '/project/mugi/nas/PAPER2/CCLM-DCEP-Tree/simulations/STD/laf2015062200_REFnew.nc'
bdc_path = '/project/mugi/nas/PAPER2/CCLM-DCEP-Tree/ahf/lbff/'

# Loading COSMO grid
nc = Dataset(cosmo_path,'r')
lon_m = nc.variables['lon'][:]
lat_m = nc.variables['lat'][:]

# Inizialization
x = np.shape(lon_m)[1]
y = np.shape(lon_m)[0]
VAR_r = np.zeros([y, x])

# list COSMO input files
f_bdc = glob.glob(bdc_path + "*.nc")
f_bdc.sort()

t = 1 # LT = UCT + 1
nday = 1
for j in range(0, ndays*24):
    data_name = name + '_' + str(t) + '_' + str(nday) + \
            '_' + str(jday) + '_' + str(sp) + '_' + str(frv) + str(ext)
    path = data_path + data_name
    # Read dataset
    print('Reading ' + data_name)
    data = gdal.Open(path)
    # Get lan and lon arrays
    lon_res = data.GetGeoTransform()[1]
    lon_l = data.GetGeoTransform()[0] + lon_res/2.1
    lon_r = data.GetGeoTransform()[0] + data.RasterXSize * lon_res - lon_res/2.1
    lon = np.arange(lon_l, lon_r, lon_res)
    lat_res = data.GetGeoTransform()[5]
    lat_l = data.GetGeoTransform()[3] + lat_res/2.1
    lat_r = data.GetGeoTransform()[3] + data.RasterYSize * lat_res - lat_res/2.1
    lat = np.arange(lat_l, lat_r, lat_res)
    lon_s = np.zeros((len(lat),len(lon)))
    for k in range(0, len(lat)):  # maybe this can be done in 1 line
        lon_s[k,:] = lon

    lat_s = np.zeros((len(lat),len(lon)))
    for k in range(0, len(lon)):  # maybe this can be done in 1 line
        lat_s[:,k] = lat

    # Read as array
    data = data.ReadAsArray()
    data = data.astype(float)
    # Prepare the data for the reshaping procedure
    lon_s_f = lon_s.flatten()
    lat_s_f = lat_s.flatten()
    data_s_f = data.flatten()
    lon_m_f = lon_m.flatten()
    lat_m_f = lat_m.flatten()
    grid_s = [lon_s_f, lat_s_f]
    grid_s = np.transpose(grid_s)
    grid_m = [lon_m_f, lat_m_f]
    grid_m = np.transpose(grid_m)
    VAR_s = griddata(grid_s, data_s_f, grid_m, method='nearest')
    VAR_s_r = np.reshape(VAR_s, (y, x))
    print('Max = ' + str(np.nanmax(VAR_s_r)))
    # Appending VAR to NetCDF file
    f = Dataset(f_bdc[j], 'a', format='NETCDF4')
    var = f.createVariable('AHF','f4',('rlat','rlon'))
    var[:] = VAR_s_r
    var.units = 'W m-2'
    var.standard_name = 'AHF'
    var.long_name = 'Anthropogenic Heat Flux'
    var.coordinates = 'lon lat'
    var.grid_mapping = 'rotated pole'
    f.close()
    # Update indexes
    if t==23:
        nday = nday + 1
        jday = jday + 1
        t = 0
    else:
        t = t + 1



