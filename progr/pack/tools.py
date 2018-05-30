#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

import numpy as np

def mask(VAR, mask):
    VAR[mask != 1] = np.nan
    VAR = np.ma.masked_invalid(VAR)
    VAR.fill_value = -999
   
    return VAR


