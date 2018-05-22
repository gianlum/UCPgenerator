#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

####################################################################
# Calculate street and building width according to Martilli 2009    
####################################################################
# INPUTS
# lp = planar area density
# lf = frontal area density
# h = mean building height
# OUTPUT
# b = building width
# w = street width
####################################################################

import numpy as np

# Constant 

def bandw(lp,lf,h):
    # Calculate building width (5a)
    b = lp/lf*h
    # Calculate street width (5b)
    w = (1/lp-1)*lp/lf*h
    print('b = ' + str(b))
    print('w = ' + str(w))
    return b, w

