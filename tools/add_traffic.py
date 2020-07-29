#!python

import glob
import numpy as np
from netCDF4 import Dataset

f_name = 'COSMOtransport_full_electric.csv'
path_f = '/net/ch4/landclim/mussettg/Singapore/d3/int2lm/lbff_traff_ev/'
path_d = '/net/ch4/landclim/mussettg/Singapore/data/traffic/'

files = glob.glob(path_f + "*.nc")
files.sort()

# extract traffic emissions
d = np.genfromtxt(path_d + f_name, delimiter=',', dtype=np.float64) # load csv file
d = d[:,2::] # discard coordinates
d = d / ((0.004 * 111000) ** 2) # convert to W m-2 (0.004 is the res in deg)
d = np.reshape(d, (200,170,24), order='C').transpose() # reshape to time, lat, lon

t = 8 # UTC + 8
for fic in files:
    print("Working on {}".format(fic))
    f = Dataset(fic, 'a')
    hf_traff = f.createVariable('HF_TRAFF','f4',('time','rlat','rlon'),
            fill_value=-1e+20)
    hf_traff[0,:,:] = d[t,:,:]
    if t == 23: # hourly loop
        t = 0
    else:
        t = t + 1
    
    # write attributed
    hf_traff.units = 'W m-2'
    hf_traff.standard_name = 'Heat flux from Traffic'
    hf_traff.long_name = 'Sensible heat flux from traffic'
    hf_traff.coordinates = 'lon lat'
    hf_traff.grid_mapping = 'rotated_pole'
    # close NetCDF
    f.close()


