#!python

######################################
# Urban Pre-processor for COSMO-DCEP #
######################################

import copy
import numpy as np
from netCDF4 import Dataset
import shapefile
import gdal

import pack.urbfrac as urbfrac
import pack.cangeom as cangeom
import pack.veget as veget
import pack.tools as tools

## USER DEFINED PATHS AND PARAMETERS

#Flags
LAD_flag = 1 # calculate LAD
threshold = 1 # for urban fraction
thr_val = 0.4 
greening = True # use greening scenario

# Height cluster
new_approach = 0 # if 1, kees fr_roof 0 at the ground
norm_vert = True # normalize roof fraction over the vertical direction

# Path to datasets
nc_path = '/project/mugi/nas/PAPER2/CCLM-DCEP-Tree/int2lm/laf2015062200.nc'
sf_path = '/project/mugi/nas/PAPER2/datasets/buildings/geom/3dbuildings_masked_geom.shp'
ufrac_path = '/project/mugi/nas/PAPER2/datasets/land_use/mosaic_20m_sealing_v2_WGS_cutted.tif'
veg_path = '/project/mugi/nas/PAPER2/datasets/trees/VEG_WGS84.tif'


## CODE (MODIFICATIONS CAN PRODUCE ERRONEOUS RESULTS)

print('Initialization')

# Open netCDF file
nc = Dataset(nc_path, 'a')

# Importing dimensions
rlon_d = len(nc.dimensions['rlon'])
rlat_d = len(nc.dimensions['rlat'])
time_d = len(nc.dimensions['time'])
udir_d = 4 
uheight1_d = 13
uclass_d = 1

# Importing variables
rlon_v = nc.variables['rlon'][:]
rlat_v = nc.variables['rlat'][:]
lon_v = nc.variables['lon'][:]
lat_v = nc.variables['lat'][:]
udir_v = np.array([-45, 0, 45, 90], dtype=np.float32)
#uheight1_v = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100], dtype=np.float32)
uheight1_v = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60], dtype=np.float32)

# Writing dimensions to the new netCDF
udir = nc.createDimension('udir',udir_d)
uheight1 = nc.createDimension('uheight1',uheight1_d)
uclass = nc.createDimension('uclass',uclass_d)

# Writing variables 
udir = nc.createVariable('udir','f4',('udir',))
uheight1 = nc.createVariable('uheight1','f4',('uheight1'))
uclass = nc.createVariable('uclass','f4',('uclass',))

udir[:] = udir_v
uheight1[:] = uheight1_v

# Initializing the variables
FR_ROOF = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
FR_URBANCL = np.zeros((1, rlat_d, rlon_d))
FR_URBAN = np.zeros((1, rlat_d, rlon_d))
M = np.zeros((1, rlat_d, rlon_d))
FR_STREETD = np.zeros((1, udir_d, rlat_d, rlon_d))
LAD_C = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
LAD_B = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
OMEGA = np.zeros((1, rlat_d, rlon_d))
LAI_URB = np.zeros((1, rlat_d, rlon_d))
mask3D = np.zeros((1, udir_d, rlat_d, rlon_d))
mask4D = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))

LAI_2 = nc.variables['LAI'][:]
Z0_2 = nc.variables['Z0'][:]
PLCOV_2 = nc.variables['PLCOV'][:]   

## CALCULATIONS

# Extracting the urban fraction
print('Extracting the urban fraction')
[FR_URBAN, data1, lat_s_1, lon_s_1] = urbfrac.wgs(ufrac_path, threshold, thr_val, rlat_v, rlon_v, FR_URBAN) 

# Calculating the street canyon geometry parameters
print('Calculating canyon parameters')
[BUILD_W, STREET_W, FR_ROOF, FR_STREETD, shapes] = cangeom.shp(sf_path, rlat_d, rlon_d, udir_d, \
        uheight1_d, rlat_v, rlon_v, FR_URBAN, FR_ROOF, norm_vert)

# Calculating the leaf area density
if LAD_flag==1:
    print('Calculating vegetation density')
    [LAD_C,OMEGA,LAI_URB] = veget.lidar(veg_path, thr_val, data1, rlat_v, rlon_v, lat_s_1, lon_s_1, LAD_C, OMEGA, LAI_URB, greening)

## MASKING NON-URBAN CELLS
print('Masking the urban cells')
# Defining the urban class mask
FR_URBANCL = copy.deepcopy(FR_URBAN)
FR_URBANCL[FR_URBANCL >= 0.11] = 1
M[0, 2:-2, 2:-2] = 1 # Mask for the cells near the boundary
FR_URBANCL = FR_URBANCL * M
FR_URBANCL[FR_URBANCL != 1] = 0
#FR_URBANCL = np.ma.masked_invalid(FR_URBANCL)
#FR_URBANCL.fill_value = -999

# Updating the mask with FR_ROOF values
FR_ROOF_mask = copy.deepcopy(np.sum(FR_ROOF,(1,2)))
FR_ROOF_mask[FR_ROOF_mask>0.01] = 1
FR_URBANCL[FR_ROOF_mask != 1] = 0

# Debug: removing grid cells with issues
FR_URBANCL[0,107-1,151-1] = 0
FR_URBANCL[0,120-1,179-1] = 0

# Updating the mask with LAD values
if LAD_flag==1:
    LAD_mask = copy.deepcopy(LAD_C[:,0,1,:,:])
    LAD_mask[LAD_mask>0.0001] = 1
    FR_URBANCL[LAD_mask != 1] = 0

mask2D = FR_URBANCL
mask3D[:,:,:,:] = FR_URBANCL[0,np.newaxis,:,:]
mask4D[:,:,:,:,:] = FR_URBANCL[0,np.newaxis,np.newaxis,:,:]

FR_URBANCL = tools.mask(FR_URBANCL, mask2D)
FR_URBAN = tools.mask(FR_URBAN, mask2D)
FR_ROOF = tools.mask(FR_ROOF, mask4D)
FR_STREETD = tools.mask(FR_STREETD, mask3D)
STREET_W = tools.mask(STREET_W, mask3D)
BUILD_W = tools.mask(BUILD_W, mask3D)
LAD_C = tools.mask(LAD_C, mask4D)
LAD_B = tools.mask(LAD_B, mask4D)
OMEGA = tools.mask(OMEGA, mask2D)

# Creating the additional "_2" variables for DCEP
PLCOV_2[FR_URBANCL == 1] = 0.8
LAI_2[FR_URBANCL == 1] = 1 # to account for grass
Z0_2[FR_URBANCL == 1] = 0.1

# Adding non-street treesa
LAI_URB[mask2D != 1] = 0 # only vegetation on urban cells
LAI_2 = LAI_2 + LAI_URB

## WRITING BACK TO NETCDF
print('Writing variables to netCDF')
# Defining the Variables
fr_roof = nc.createVariable('FR_ROOF','f4',('uclass','udir','uheight1','rlat','rlon',),fill_value=-1e+20)
fr_urbancl = nc.createVariable('FR_UCLASS','f4',('uclass','rlat','rlon'),fill_value=-1e+20)
fr_urban = nc.createVariable('FR_URB','f4',('uclass','rlat','rlon'),fill_value=-1e+20)
fr_streetd = nc.createVariable('FR_UDIR','f4',('uclass','udir','rlat','rlon',),fill_value=-1e+20)
street_w = nc.createVariable('W_STREET','f4',('uclass','udir','rlat','rlon',),fill_value=-1e+20)
build_w = nc.createVariable('W_BUILD','f4',('uclass','udir','rlat','rlon'),fill_value=-1e+20)
lad_c = nc.createVariable('LAD_C','f4',('uclass','udir','uheight1','rlat','rlon'),fill_value=-1e+20)
lad_b = nc.createVariable('LAD_B','f4',('uclass','udir','uheight1','rlat','rlon'),fill_value=-1e+20)
omega_r = nc.createVariable('OMEGA_R','f4',('uclass','rlat','rlon'),fill_value=-1e+20)
omega_d = nc.createVariable('OMEGA_D','f4',('uclass','rlat','rlon'),fill_value=-1e+20)

lai_2 = nc.createVariable('LAI_2','f4',('time','rlat','rlon'),fill_value=-1e+20)
plcov_2 = nc.createVariable('PLCOV_2','f4',('time','rlat','rlon'),fill_value=-1e+20)
z0_2 = nc.createVariable('Z0_2','f4',('time','rlat','rlon'),fill_value=-1e+20)

# Writing attributes
fr_roof.units = '1'
fr_roof.standard_name = 'Wall surfaces fraction'
fr_roof.long_name = 'Fraction of wall surfaces per direction and height'
fr_roof.coordinates = 'lon lat'
fr_roof.grid_mapping = 'rotated_pole'
fr_roof._FillValue = -1e+20

fr_urbancl.units = '1'
fr_urbancl.standard_name = 'Urban Mask'
fr_urbancl.long_name = 'Mask for urban area'
fr_urbancl.coordinates = 'lon lat'
fr_urbancl.grid_mapping = 'rotated_pole'
fr_urbancl._FillValue = -1e+20

fr_urban.units = '1'
fr_urban.standard_name = 'Urban Fraction'
fr_urban.long_name = 'Fraction of urban horizontal surfaces'
fr_urban.coordinates = 'lon lat'
fr_urban.grid_mapping = 'rotated_pole'
fr_urban._FillValue = -1e+20

fr_streetd.units = '1'
fr_streetd.standard_name = 'Canyon Direction Fraction'
fr_streetd.long_name = 'Fraction of canyon directions'
fr_streetd.coordinates = 'lon lat'
fr_streetd.grid_mapping = 'rotated_pole'
fr_streetd._FillValue = -1e+20

street_w.units = 'm'
street_w.standard_name = 'Street Width'
street_w.long_name = 'Average street width per direction'
street_w.coordinates = 'lon lat'
street_w.grid_mapping = 'rotated_pole'
street_w._FillValue = -1e+20

build_w.units = 'm'
build_w.standard_name = 'Building Width'
build_w.long_name = 'Average building width per direction'
build_w.coordinates = 'lon lat'
build_w.grid_mapping = 'rotated_pole'
build_w._FillValue = -1e+20

lad_c.units = 'm2m-3'
lad_c.standard_name = 'Canyon LAD'
lad_c.long_name = 'Leaf Area Density in the Canyon'
lad_c.coordinates = 'lon lat'
lad_c.grid_mapping = 'rotated_pole'
lad_c._FillValue = -1e+20

lad_b.units = 'm2m-3'
lad_b.standard_name = 'Building LAD'
lad_b.long_name = 'Leaf Area Density in the Building Column'
lad_b.coordinates = 'lon lat'
lad_b.grid_mapping = 'rotated_pole'
lad_b._FillValue = -1e+20

omega_r.units = '1'
omega_r.standard_name = 'Clumping rad.'
omega_r.long_name = 'Neighborhood-scale clumping coeff. for rad. for tree foliage'
omega_r.coordinates = 'lon lat'
omega_r.grid_mapping = 'rotated_pole'
omega_r._FillValue = -1e+20

omega_d.units = '1'
omega_d.standard_name = 'Clumping drag'
omega_d.long_name = 'Neighborhood-scale clumping coeff. for drag for tree foliage'
omega_d.coordinates = 'lon lat'
omega_d.grid_mapping = 'rotated_pole'
omega_d._FillValue = -1e+20

lai_2.units = 'm2m-2'
lai_2.standard_name = 'Leaf Area Index 2'
lai_2.long_name = 'Leaf Area Index for urban vegetation'
lai_2.coordinates = 'lon lat'
lai_2.grid_mapping = 'rotated_pole'
lai_2._FillValue = -1e+20

plcov_2.units = '1'
plcov_2.standard_name = 'Plant Coverage 2'
plcov_2.long_name = 'Plant coverage for urban vegetation'
plcov_2.coordinates = 'lon lat'
plcov_2.grid_mapping = 'rotated_pole'
plcov_2._FillValue = -1e+20

z0_2.units = 'm'
z0_2.standard_name = 'Roughness Length 2' 
z0_2.long_name = 'Roughness Length for urban vegetation'
z0_2.coordinates = 'lon lat'
z0_2.grid_mapping = 'rotated_pole'
z0_2._FillValue = -1e+20

# Inserting data into variables
fr_roof[:] = FR_ROOF
fr_urbancl[:] = FR_URBANCL
fr_urban[:] = FR_URBAN
fr_streetd[:] = FR_STREETD
street_w[:] = STREET_W
build_w[:] = BUILD_W
lad_c[:] = LAD_C
lad_b[:] = LAD_B
omega_r[:] = OMEGA
omega_d[:] = OMEGA

lai_2[:] = LAI_2
plcov_2[:] = PLCOV_2
z0_2[:] = Z0_2

# Closing the netCDF
nc.close()

print ('New netCDF file generated')
