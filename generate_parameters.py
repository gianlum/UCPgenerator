#!python

######################################
# Urban Pre-processor for COSMO-DCEP #
######################################

import copy
import numpy as np
from netCDF4 import Dataset
import shapefile
import cluster
from scipy import interpolate
import gdal
import geometry
import angle
import geo2rot

#Flags
LAD_flag = 0 # calculate LAD

# Height cluster
new_approach = 0 # if 1, kees fr_roof 0 at the ground

# Opening the datasets
nc = Dataset('/project/mugi/nas/PAPER2/CCLM-DCEP-Tree/int2lm/laf2015062200.nc','a')
sf = shapefile.Reader("/project/mugi/nas/PAPER2/datasets/buildings/geom/3dbuildings_masked_geom.shp")
shapes = sf.shapes()
ufrac_path = '/project/mugi/nas/PAPER2/datasets/land_use/mosaic_20m_sealing_v2_WGS_cutted.tif'
veg_path = '/project/mugi/nas/PAPER2/datasets/trees/VEG_WGS84.tif'

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
AREA_BLD = np.zeros((1, rlat_d, rlon_d))
VERT_AREA = np.zeros((1, udir_d, rlat_d, rlon_d))
MEAN_HEIGHT = np.zeros((1, udir_d, rlat_d, rlon_d))
FR_ROOF = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
FR_URBANCL = np.zeros((1, rlat_d, rlon_d))
FR_URBAN = np.zeros((1, rlat_d, rlon_d))
FR_URBAN_count = np.zeros((1, rlat_d, rlon_d))
M = np.zeros((1, rlat_d, rlon_d))
FR_STREETD = np.zeros((1, udir_d, rlat_d, rlon_d))
FR_STREET = np.zeros((1, udir_d, rlat_d, rlon_d))
STREET_W = np.zeros((1, udir_d, rlat_d, rlon_d))
BUILD_W = np.zeros((1, udir_d, rlat_d, rlon_d))
LAD_C = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
LAD_B = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))
OMEGA = np.zeros((1, rlat_d, rlon_d))
LAI_URB = np.zeros((1, rlat_d, rlon_d))
mask = np.zeros((1, udir_d, rlat_d, rlon_d))
mask2 = np.zeros((1, udir_d, uheight1_d, rlat_d, rlon_d))

LAI_2 = nc.variables['LAI'][:]
Z0_2 = nc.variables['Z0'][:]
PLCOV_2 = nc.variables['PLCOV'][:]   

# Extracting the urban fraction
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
        if data1_tmp >= 0.5:
            FR_URBAN[0,lat_idx,lon_idx] += 1.
    


FR_URBAN = FR_URBAN/FR_URBAN_count

# Defining the urban class mask
FR_URBANCL = copy.deepcopy(FR_URBAN)
FR_URBANCL[FR_URBANCL >= 0.11] = 1
M[0, 2:-2, 2:-2] = 1 # Mask for the cells near the boundary
FR_URBANCL = FR_URBANCL * M
FR_URBANCL[FR_URBANCL != 1] = 0
#FR_URBANCL = np.ma.masked_invalid(FR_URBANCL)
#FR_URBANCL.fill_value = -999

# Calculating the building area (vertical and horizontal)
N = len(shapes)
for x in range (0, N):
    # Reading Geometry
    shape = shapes[x]
    p = shape.points
    p = np.array(p)
    # Reading the Area (horizontal)
    area = sf.record(x)[4]
    # Calculating the coordinates of the centroid of the polygon
    lonC = np.mean(p[:,0])
    latC = np.mean(p[:,1])
    # Converting centroid to rotated coordinates
    [lonC_r,latC_r] = geo2rot.g2r(lonC,latC)
    # Calculating the index of the correspoding grid point
    lon_idx = np.abs(rlon_v - lonC_r).argmin()
    lat_idx = np.abs(rlat_v - latC_r).argmin()
    # Allocating the total building area
    AREA_BLD[0,lat_idx,lon_idx]+=area
    # Clustering the geomerty heights
    hgt = sf.record(x)[0]
    MEAN_HEIGHT[0,:,lat_idx,lon_idx] += hgt * area
    if new_approach==1:
        hgt_class = cluster.height_new(hgt)
    else:
        hgt_class = cluster.height_old(hgt)

    for k in range (1, len(p)): # Looping over building segments (facades)
        vert_area = geometry.distWGS(p,k) * hgt 
        ang = geometry.angle(p,k)
        ang_class = cluster.angle(ang)
        # 
        VERT_AREA[0,ang_class,lat_idx,lon_idx] += vert_area
        # Allocating the building area by direction and height
        FR_ROOF[0,ang_class,hgt_class,lat_idx,lon_idx]+=vert_area


# Calculating the area densities (Grimmond and Oke, 1998)
area_grid = abs(rlon_v[1]-rlon_v[0])*abs(rlat_v[1]-rlat_v[0]) # regular grida
lambda_p = AREA_BLD/area_grid*FR_URBAN
lambda_f = VERT_AREA/area_grid*FR_URBAN

# Calculating the street and building widths (Martilli, 2009)
h_m = MEAN_HEIGHT / AREA_BLD
BUILD_W = lambda_p / lambda_f * h_m
STREET_W = (1 / lambda_p - 1) * lambda_p / lambda_f * h_m

if LAD_flag==1:
    print('Calculating LAD')
    # Read the dataset
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
    # Write LAD into cosmo output
    for j in range(0, len(lon)):
        for k in range(0, len(lat)):
            data2_tmp = data2[k,j]
            if data2_tmp > 0.:
                # Identify corresponding grid cell
                [lonV,latV] = np.array(geo2rot.g2r(lon_s_2[k,j],lat_s_2[k,j]))
                lon_idx = np.abs(rlon_v - lonC).argmin()
                lat_idx = np.abs(rlat_v - latC).argmin()
                # Calculate LAD from Lidar canopy height
                area_lidar = 2.23 # m2 mesured in GIS
                LAD_spec = 1. # m2 m-3 (Klingberg et al., 2017)
                area_grid = 278. * 278. # m2
                h_ug = 5. # m, height urban grida
                h_cbase = 3. # m, height of the canopy base
                data2_tmp = data2[k,j]
                # Check if in-canyon vegetation
                lon_urb_idx = np.abs(lon_s_1[0,:] - lon_s_2[k,j]).argmin()
                lat_urb_idx = np.abs(lat_s_1[:,0] - lat_s_2[k,j]).argmin()
                data1_tmp = data1[lat_urb_idx,lon_urb_idx]
                if data1_tmp >= 0.5:
                    # In-canyon vegetation
                    if data2_tmp >= h_cbase : 
                        LAD_C[0,:,0,lat_idx,lon_idx] += LAD_spec * area_lidar * \
                            (h_ug - h_cbase) / area_grid / h_ug
                    if data2_tmp >= 5 :
                        LAD_C[0,:,1,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid 
                    if data2_tmp >= 10 :
                        LAD_C[0,:,2,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                    if data2_tmp >= 15 :
                        LAD_C[0,:,3,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                    if data2_tmp >= 20 :
                        LAD_C[0,:,4,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                    if data2_tmp >= 25 :
                        LAD_C[0,:,5,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                    if data2_tmp >= 30 :
                        LAD_C[0,:,6,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                    if data2_tmp >= 35 :
                        LAD_C[0,:,7,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                    if data2_tmp >= 45 :
                        LAD_C[0,:,8,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                    if data2_tmp >= 50 :
                        LAD_C[0,:,9,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                    if data2_tmp >= 55 :
                        LAD_C[0,:,10,lat_idx,lon_idx] += LAD_spec * area_lidar / area_grid
                else:
                    # Out-canyon vegetation
                    LAI_URB[0,lat_idx,lon_idx] += LAD_spec * area_lidar * data2_tmp / area_grid






# Calculating the foliage clumping coefficient (OMEGA)
# TO DO: implement a way to calculate from tree distribution
OMEGA[0,:,:] = 0.5 # more relistic value for now

# Calculating and normalizing the canyon direction distribution
FR_STREETD = np.sum(FR_ROOF,2)
norm_streetd = np.sum(FR_STREETD,1)
FR_STREETD = FR_STREETD/norm_streetd

# Normalizing FR_ROOF
norm_fr_roof = np.sum(FR_ROOF,1)
for j in range(0, udir_d):
    for k in range(0, uheight1_d):
        for o in range(0, rlat_d):
            for z in range(0, rlon_d):
                if norm_fr_roof[0,k,o,z] != 0:
                    FR_ROOF[0,j,k,o,z] = FR_ROOF[0,j,k,o,z]/norm_fr_roof[0,k,o,z]





FR_ROOF[FR_ROOF<0.001] = 0 # to avoid negative values 

# Updating the mask with LAD values

if LAD_flag:
    LAD_mask = copy.deepcopy(LAD_C[:,0,1,:,:])
    LAD_mask[LAD_mask>0] = 1
    FR_URBANCL[LAD_mask != 1] = 0

# Masking the 0 values # TO COMPLETE AT THE END
mask[:,:,:,:] = FR_URBANCL[0,np.newaxis,:,:]
mask2[:,:,:,:,:] = FR_URBANCL[0,np.newaxis,np.newaxis,:,:]

FR_URBAN[FR_URBANCL != 1] = np.nan
FR_URBAN = np.ma.masked_invalid(FR_URBAN)
FR_URBAN.fill_value = -999

FR_ROOF[mask2 != 1] = np.nan
FR_ROOF = np.ma.masked_invalid(FR_ROOF)
FR_ROOF.fill_value = -999

FR_STREETD[mask != 1] = np.nan
FR_STREETD = np.ma.masked_invalid(FR_STREETD)
FR_STREETD.fill_value = -999

STREET_W[mask != 1] = np.nan
STREET_W = np.ma.masked_invalid(STREET_W)
STREET_W.fill_value = -999

BUILD_W[mask != 1] = np.nan
BUILD_W = np.ma.masked_invalid(BUILD_W)
BUILD_W.fill_value = -999

LAD_C[mask2 != 1] = np.nan
LAD_C = np.ma.masked_invalid(LAD_C)
LAD_C.fill_value = -999

LAD_B[mask2 != 1] = np.nan
LAD_B = np.ma.masked_invalid(LAD_B)
LAD_B.fill_value = -999

OMEGA[FR_URBANCL != 1] = np.nan
OMEGA = np.ma.masked_invalid(OMEGA)
OMEGA.fill_value = -999

# Creating the additional "_2" variables for DCEP
PLCOV_2[FR_URBANCL == 1] = 0.8
LAI_2[FR_URBANCL == 1] = 3
Z0_2[FR_URBANCL == 1] = 0.1

# Defining the Variables
fr_roof = nc.createVariable('FR_ROOF','f4',('uclass','udir','uheight1','rlat','rlon',),fill_value=-999)
fr_urbancl = nc.createVariable('FR_UCLASS','f4',('uclass','rlat','rlon'),fill_value=-999)
fr_urban = nc.createVariable('FR_URB','f4',('uclass','rlat','rlon'),fill_value=-999)
fr_streetd = nc.createVariable('FR_UDIR','f4',('uclass','udir','rlat','rlon',),fill_value=-999)
street_w = nc.createVariable('W_STREET','f4',('uclass','udir','rlat','rlon',),fill_value=-999)
build_w = nc.createVariable('W_BUILD','f4',('uclass','udir','rlat','rlon'),fill_value=-999)
lad_c = nc.createVariable('LAD_C','f4',('uclass','udir','uheight1','rlat','rlon'),fill_value=-999)
lad_b = nc.createVariable('LAD_B','f4',('uclass','udir','uheight1','rlat','rlon'),fill_value=-999)
omega_r = nc.createVariable('OMEGA_R','f4',('uclass','rlat','rlon'),fill_value=-999)
omega_d = nc.createVariable('OMEGA_D','f4',('uclass','rlat','rlon'),fill_value=-999)
lai_urban = nc.createVariable('LAI_URB','f4',('uclass','rlat','rlon'),fill_value=-999)

lai_2 = nc.createVariable('LAI_2','f4',('time','rlat','rlon'))
plcov_2 = nc.createVariable('PLCOV_2','f4',('time','rlat','rlon'))
z0_2 = nc.createVariable('Z0_2','f4',('time','rlat','rlon'))

# Writing attributes
fr_roof.units = '1'
fr_roof.standard_name = 'Wall surfaces fraction'
fr_roof.long_name = 'Fraction of wall surfaces per direction and height'
fr_roof.coordinates = 'lon lat'
fr_roof.grid_mapping = 'rotated_pole'

fr_urbancl.units = '1'
fr_urbancl.standard_name = 'Urban Mask'
fr_urbancl.long_name = 'Mask for urban area'
fr_urbancl.coordinates = 'lon lat'
fr_urbancl.grid_mapping = 'rotated_pole'

fr_urban.units = '1'
fr_urban.standard_name = 'Urban Fraction'
fr_urban.long_name = 'Fraction of urban horizontal surfaces'
fr_urban.coordinates = 'lon lat'
fr_urban.grid_mapping = 'rotated_pole'

fr_streetd.units = '1'
fr_streetd.standard_name = 'Canyon Direction Fraction'
fr_streetd.long_name = 'Fraction of canyon directions'
fr_streetd.coordinates = 'lon lat'
fr_streetd.grid_mapping = 'rotated_pole'

street_w.units = 'm'
street_w.standard_name = 'Street Width'
street_w.long_name = 'Average street width per direction'
street_w.coordinates = 'lon lat'
street_w.grid_mapping = 'rotated_pole'

build_w.units = 'm'
build_w.standard_name = 'Building Width'
build_w.long_name = 'Average building width per direction'
build_w.coordinates = 'lon lat'
build_w.grid_mapping = 'rotated_pole'

lad_c.units = 'm2m-3'
lad_c.standard_name = 'Canyon LAD'
lad_c.long_name = 'Leaf Area Density in the Canyon'
lad_c.coordinates = 'lon lat'
lad_c.grid_mapping = 'rotated_pole'

lad_b.units = 'm2m-3'
lad_b.standard_name = 'Building LAD'
lad_b.long_name = 'Leaf Area Density in the Building Column'
lad_b.coordinates = 'lon lat'
lad_b.grid_mapping = 'rotated_pole'

omega_r.units = '1'
omega_r.standard_name = 'Clumping rad.'
omega_r.long_name = 'Neighborhood-scale clumping coeff. for rad. for tree foliage'
omega_r.coordinates = 'lon lat'
omega_r.grid_mapping = 'rotated_pole'

omega_d.units = '1'
omega_d.standard_name = 'Clumping drag'
omega_d.long_name = 'Neighborhood-scale clumping coeff. for drag for tree foliage'
omega_d.coordinates = 'lon lat'
omega_d.grid_mapping = 'rotated_pole'

lai_urban.units = '1'
lai_urban.standard_name = 'Out-Canyon Leaf Area Index'
lai_urban.long_name = 'Leaf Area Index from vegetation outside the street canyon'
lai_urban.coordinates = 'lon lat'
lai_urban.grid_mapping = 'rotated_pole'

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

z0_2.units = 'm'
z0_2.standard_name = 'Roughness Length 2' 
z0_2.long_name = 'Roughness Length for urban vegetation'
z0_2.coordinates = 'lon lat'
z0_2.grid_mapping = 'rotated_pole'

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
lai_urban[:] = LAI_URB

lai_2[:] = LAI_2
plcov_2[:] = PLCOV_2
z0_2[:] = Z0_2

# Closing the netCDF
nc.close()

print ('New netCDF file generated')
