
DEBUG=True

import copy
import numpy as np
from netCDF4 import Dataset
import shapefile
import find
import cluster
from scipy import interpolate

if DEBUG:
    import matplotlib.pyplot as plt


import geometry
import angle
import geo2rot


##
nc = Dataset('/media/pick/Data/DCEP/laf2015062200_intermediate.nc','a')
nc_lu = Dataset('/media/pick/Data/DCEP/basel_corine_reproject.nc','r')



rlon_d = len(nc.dimensions['rlon'])
rlat_d = len(nc.dimensions['rlat'])
time_d = len(nc.dimensions['time'])
udir_d = 4 
uheight1_d = 13
uclass_d = 1

rlon_v = nc.variables['rlon'][:]
rlat_v = nc.variables['rlat'][:]
lon_v = nc.variables['lon'][:]
lat_v = nc.variables['lat'][:]
FR_URBAN = nc.variables['FR_URB'][:]
lon_v_lu = nc_lu.variables['lon'][:]
lat_v_lu = nc_lu.variables['lat'][:]
var_lu = nc_lu.variables['Band1'][:]
udir_v = np.array([-45, 0, 45, 90], dtype=np.float32)
#uheight1_v = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100], dtype=np.float32)
uheight1_v = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], dtype=np.float32)

AREA_BLD = np.zeros((1, rlat_d, rlon_d))
FR_ROOF = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
FR_ROOF2 = np.zeros((1, uheight1_d, rlat_d, rlon_d))
FR_URBANCL = np.zeros((1, rlat_d, rlon_d))
M = np.zeros((1, rlat_d, rlon_d))
FR_STREETD = np.zeros((1, udir_d, rlat_d, rlon_d))
FR_STREET = np.zeros((1, udir_d, rlat_d, rlon_d))
STREET_W = np.zeros((1, udir_d, rlat_d, rlon_d))
BUILD_W = np.zeros((1, udir_d, rlat_d, rlon_d))
LAD_C = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
LAD_B = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
OMEGA = np.zeros((1, rlat_d, rlon_d))
mask = np.zeros((1, udir_d, rlat_d, rlon_d))
mask2 = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))


xx, yy = np.meshgrid(lon_v_lu, lat_v_lu)

coords_lu = np.append(xx.flatten()[:,np.newaxis], yy.flatten()[:,np.newaxis], axis = 1)
var_lu = var_lu.flatten()
# Meshgrid for rotated latlon coordinates
RLON_V, RLAT_V = np.meshgrid(rlon_v, rlat_v)
coords_mod = np.append(RLON_V.flatten()[:,np.newaxis], RLAT_V.flatten()[:,np.newaxis], axis = 1)

VAR = interpolate.griddata(coords_lu, var_lu, coords_mod, method = 'nearest')

if DEBUG:
    print(FR_URBAN[0])
    plt.contour(FR_URBAN[0])
    plt.show()
    area_grid = abs(rlon_v[1]-rlon_v[0])*abs(rlat_v[1]-rlat_v[0])
    print(area_grid)