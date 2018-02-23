#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

######################################
# 14.10.2016 Gianluca Mussetti
# Read the shapefile in WGS 84, rotate it, calculate the building orientation and
# write down in a new shapefile
# Written combining poly_shp_mod.py, poly_shp_dir.py and adapting to the 
# swissBUILDINGS3D_1.0 (LoD1 building geometry dataset for all Switzerland)
######################################

# BUILT-IN MODULES
import sys
import os
import string
import pickle
import datetime
import struct
import numpy as np
from netCDF4 import Dataset
# OTHER MODULES
import shapefile
import geo2rot
import angle
#import length
import perimeter

###########
# READING # 
###########

# Reading the shapefile
sf = shapefile.Reader("/home/mugi/DATA/WGS84/Buildings/swissBUILDINGS3D_1.0/buildings3D.shp")
shapes = sf.shapes()

# Open a new shapefile
w = shapefile.Writer(shapefile.POLYGON)

# Reading the Fields (constants)
fields = sf.fields
field1 = fields[6] # Building Height

# Writing the Fields
w.field(field1[0],field1[1],field1[2],field1[3])
#w.field('Canyon Dir','N',5,1)
w.field('Vert_Area','N',47,15) # Building Surface (vertical)
# not sure about this
w.field('Horiz_Area','N',47,15) # Building Surface (horizontal)

# Iterating over N
N = len(shapes)
for x in range (0, N):
   # Reaging Geometry
   shape = shapes[x]
   p = shape.points
   p = np.array(p)
   lon_g = p[:,0]
   lat_g = p[:,1]
   # Reading the Records
   records = sf.record(x)
   rec1 = records[-1]
#   # Calculating the angle
#   seg = length.dist(p)
#   ang = angle.alpha(seg)
#   rec2 = round(ang,1)
   # Caluculating the vertical surface area
   perim = perimeter.per(p)
   area = perim * rec1
   # Tranforming the coordinates
   [lon_r,lat_r] = geo2rot.g2r(lon_g,lat_g)
   # Calculating the Area
   rec3 = 0.5*np.abs(np.dot(lon_r,np.roll(lat_r,1))-np.dot(lat_r,np.roll(lon_r,1)))
   # Caluculating the vertical surface area
   perim = perimeter.per(p)
   rec2 = perim * rec1
   # Defining back the Geometry
   p_r = [lon_r,lat_r]
   p_r = np.transpose(p_r)
   p_r = np.array(p_r).tolist()
   # Writing the Geometry
   w.poly(parts=[p_r])
   w.record(rec1,rec2,rec3)

# Saving the new shapefile
w.save('/home/mugi/DATA/ROTATED/Buildings/swissBUILDINGS3D_1.0/buildings3D')


