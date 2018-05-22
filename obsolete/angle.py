#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

##################################################################################
# Calculate the angle of the street canyon, as defined 
# in my notes. Basically is the angle clockwise from the north 
# It has values from 0 to 180 degree.
# INPUT: numpy array containing the coordinates of the longest vector
# OUTPUT: numpy array containing the (1) angle as defined above
##################################################################################

import numpy as np
#import find

def alpha(array):
    y1 = array[0,1]
    y2 = array[1,1]
    x1 = array[0,0]
    x2 = array[1,0]
    x = array[:,0]
    alpha_f = np.degrees(np.arctan2([y2-y1],[x2-x1]))
    # identification of the case
    if (x2>x1 and y2>y1):
        alpha = 90 - alpha_f
    elif (x2>x1 and y2<y1):
        alpha = 90 + alpha_f
    elif (x2<x1 and y2>y1):
        alpha = 270 - alpha_f
    else:
        alpha = abs(alpha_f) - 90

    return alpha

def categ(angle):
    # identification of the category
    if (angle<22.5):
        case = 'N-S'
    elif (angle<67.5):
        case = 'NE-SW'
    elif (angle<112.5):
        case = 'W-E'
    elif (angle<157.5):
        case = 'NW-SE'
    else:
        case = 'N-S'

    return case
