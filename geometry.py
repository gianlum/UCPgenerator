#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

#################################################################
# Module that includes some basics model calculations           #
# dist -> calculates the distance between point k and point k+1 #
# 
#################################################################

import numpy as np

def dist(array,k):
    x1 = array[k-1][0]
    x2 = array[k][0]
    y1 = array[k-1][1]
    y2 = array[k][1]
    dist = ((x2-x1)**2 + (y2-y1)**2)**0.5 # not in m
    return dist

def distWGS(array,k):
    # Using the haversine formula
    R = 6371 # radius of earth in km
    x1 = array[k-1][0]
    x2 = array[k][0]
    y1 = array[k-1][1]
    y2 = array[k][1]a
    dLon = np.deg2rad(x2-x1)
    dLat = np.deg2rad(y2-y1)
    a = np.sin(dLat/2.) * np.sin(dLat/2.) + \
            np.cos(np.deg2rad(y1)) * np.cos(np.deg2rad(y2)) * \
            np.sin(dLon/2.) * np.sin(dLon)
    c = 2. * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    dist = R * c * 1000 # distance in m
    return dist

def angle(array,k):
    y1 = array[k-1][1]
    y2 = array[k][1]
    x1 = array[k-1][0]
    x2 = array[k][0]
    alpha_f = np.degrees(np.arctan2([y2-y1],[x2-x1]))
    # identification of the case
    if (x2>=x1 and y2>=y1):
       alpha = 90 - alpha_f
    elif (x2>x1 and y2<y1):
       alpha = 90 + alpha_f
    elif (x2<x1 and y2>y1):
       alpha = 270 - alpha_f
    else:
       alpha = abs(alpha_f) - 90

    return alpha

