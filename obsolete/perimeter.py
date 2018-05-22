#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

##################################################################################
# Returns the perimeter of the polygon
# INPUT: array containing the couples of coordinated
# OUTPUT: the perimeter
##################################################################################

import numpy as np
import find

def per(array):
    array1 = array[0:len(array)-1]
    array2 = array[1:len(array)]
    dist = (np.sum((array1-array2)**2, axis=1))**0.5 # not in m
    perim = np.sum(dist)
    return perim


