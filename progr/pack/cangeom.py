#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

import numpy as np
import shapefile

from . import cluster
from . import geometry
from . import geo2rot

def shp(sf_path, rlat_d, rlon_d, udir_d, uheight1_d, rlat_v, rlon_v, \
        FR_URBAN, FR_ROOF, norm_vert):
    # Inizializations
    AREA_BLD = np.zeros((1, rlat_d, rlon_d))
    BUILD_W = np.zeros((1, udir_d, rlat_d, rlon_d))
    FR_STREETD = np.zeros((1, udir_d, rlat_d, rlon_d))
    LAMBDA_F = np.zeros((1, udir_d, rlat_d, rlon_d))
    LAMBDA_P = np.zeros((1, rlat_d, rlon_d))
    MEAN_HEIGHT = np.zeros((1, rlat_d, rlon_d))    
    STREET_W = np.zeros((1, udir_d, rlat_d, rlon_d))
    VERT_AREA = np.zeros((1, udir_d, rlat_d, rlon_d))
    # Read the dataset
    sf = shapefile.Reader(sf_path)
    shapes = sf.shapes()
    # Calculations
    N = len(shapes)
    for x in range (0, N):
        # Reading Geometry
        shape = shapes[x]
        p = shape.points
        p = np.array(p)
        # Reading the Area (horizontal)
        area = sf.record(x)[19]
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
        hgt = sf.record(x)[15]
        MEAN_HEIGHT[0,lat_idx,lon_idx] += hgt * area
        hgt_class = cluster.height_old(hgt)
        # Looping over building segments (facades) 
        for k in range (1, len(p)): 
            vert_area = geometry.distWGS(p,k) * hgt
            ang = geometry.angle(p,k)
            ang_class = cluster.angle(ang)
            # 
            VERT_AREA[0,ang_class,lat_idx,lon_idx] += vert_area
            # Allocating the building area by direction and height
            FR_ROOF[0,ang_class,hgt_class,lat_idx,lon_idx]+=vert_area
            

    # Calculating the area densities (Grimmond and Oke, 1998)
    np.seterr(divide='ignore') # disabled division by 0 warning
    area_grid = 197136 # m2, calculated in GIS. TO DO: calculate from grid
    lambda_p = AREA_BLD[0,:,:]/(area_grid*FR_URBAN[0,:,:])
    lambda_p[lambda_p>0.9] = 0.9 # upper limit
    lambda_p[lambda_p<0.1] = 0.1 # lower limit
    lambda_f = VERT_AREA[0,:,:,:]/(area_grid*FR_URBAN[0,:,:])
    lambda_f[lambda_f>0.9] = 0.9 # upper limit
    lambda_f[lambda_f<0.1] = 0.1 # lower limit 
    LAMBDA_P[0,:,:] = lambda_p
    LAMBDA_F[0,:,:,:] = lambda_f
    # Calculating the street and building widths (Martilli, 2009)
    h_m = MEAN_HEIGHT[0,:,:] / AREA_BLD[0,:,:]
    #h_m[h_m==15] = 10.
    BUILD_W[0,:,:,:] = lambda_p[np.newaxis,:,:] / lambda_f[:,:,:] * h_m
    STREET_W[0,:,:,:] = (1 / lambda_p[np.newaxis,:,:] - 1) * lambda_p[np.newaxis,:,:] \
                        / lambda_f * h_m
    #STREET_W[STREET_W<5] = 5  # min aspect ration LCZ 2
    #STREET_W[STREET_W>100] = 100  # max aspect ration LCZ 9
   
    # Calculating and normalizing the canyon direction distribution
    FR_STREETD = np.sum(FR_ROOF,2)
    norm_streetd = np.sum(FR_STREETD,1)
    FR_STREETD = FR_STREETD/norm_streetd

    # Normalizing FR_ROOF
    if norm_vert == True:
        norm_fr_roof = np.sum(FR_ROOF,2)
    else:
        norm_fr_roof = np.sum(FR_ROOF,1)
    
    for j in range(0, udir_d):
        for k in range(0, uheight1_d):
            for o in range(0, rlat_d):
                for z in range(0, rlon_d):
                    if norm_vert == True:
                        if norm_fr_roof[0,j,o,z] != 0:
                            FR_ROOF[0,j,k,o,z] = FR_ROOF[0,j,k,o,z]/norm_fr_roof[0,j,o,z]

                    else:
                        if norm_fr_roof[0,k,o,z] != 0:
                            FR_ROOF[0,j,k,o,z] = FR_ROOF[0,j,k,o,z]/norm_fr_roof[0,k,o,z]






    FR_ROOF[FR_ROOF<0.001] = 0 # to avoid negative values 

    return BUILD_W, STREET_W, FR_ROOF, FR_STREETD, shapes


