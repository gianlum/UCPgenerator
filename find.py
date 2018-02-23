#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

##################################################################################
# add the comment                             #
##################################################################################

import numpy as np

def nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
