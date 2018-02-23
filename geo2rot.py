#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

##################################################################################
# Performs the transformation from Geographical coordinate to rotated coordinate #
# lon_g = longitude array in the geographical system                             #
# lat_g = latitude array in the geographical system                              #
##################################################################################

import numpy as np

# Constant related to the user projection (from the COSMO namelist)
pollat = np.deg2rad(43.0)    # latitude of the rotated pole
pollon = np.deg2rad(-170.0)  # longitude of the rotated pole

def g2r(lon_g,lat_g):
    # Transform the coordinates in radiants
    lon_g = np.deg2rad(lon_g)
    lat_g = np.deg2rad(lat_g)
    # Apply the transformation (formula (3.74) COSMO Documentation - Dynamics and Numerics)
    lon_r = np.arctan ( ( np.cos(lat_g) * np.sin(lon_g-pollon) ) / ( np.cos(lat_g) * np.sin(pollat) * np.cos(lon_g-pollon) - np.sin(lat_g) * np.cos(pollat) ) )
    lat_r = np.arcsin( np.sin(lat_g) * np.sin(pollat) + np.cos(lat_g) * np.cos(pollat) * np.cos(lon_g-pollon) )
    # Transform new coordinates to degrees
    lon_r = np.rad2deg(lon_r)
    lat_r = np.rad2deg(lat_r)
    return lon_r, lat_r

