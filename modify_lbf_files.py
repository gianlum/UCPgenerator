#!python

import os
import numpy as np
from netCDF4 import Dataset

path = 'location_of_input_files_(lbff*)_excluded_the_laf'

files = [ fic for fic in os.listdir(path)]
files.sort()

nc_dcep = Dataset('location_of_laf_file','r')
FR_URBANCL = nc_dcep.variables['FR_UCLASS'][:]

for fic in files:
    obj = path + fic
    f = Dataset(obj, 'a', format='NETCDF4')
    lai = f.variables['LAI'][:]
    plcov = f.variables['PLCOV'][:]
    # Creating new variables
    lai_2 = f.createVariable('LAI_2','f4',('time','rlat','rlon'))
    plcov_2 = f.createVariable('PLCOV_2','f4',('time','rlat','rlon'))
    # Values of the new variables
    plcov[FR_URBANCL == 1] = 0.8
    lai[FR_URBANCL == 1] = 3
    # Inserting data into variables
    lai_2[:] = lai
    plcov_2[:] = plcov
    # Writing attributes
    lai_2.units = 'm2m-2'
    lai_2.standard_name = 'Leaf Area Index 2'
    lai_2.long_name = 'Leaf Area Index for urban vegetation'
    lai_2.coordinates = 'lon lat'
    lai_2.grid_mapping = 'rotated_pole'
    plcov_2.units = '1'
    plcov_2.standard_name = 'Plant Coverage 2'
    plcov_2.long_name = 'Plant coverage for urban vegetation'
    plcov_2.coordinates = 'lon lat'
    plcov_2.grid_mapping = 'rotated_pole'
    # Closing netCDF
    f.close()

nc_dcep.close()

print 'lbff files modified'
