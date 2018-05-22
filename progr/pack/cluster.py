#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

##################################################################################
##################################################################################

import numpy as np
#import find
# here you should than start the loop

def angle(angle):
    # identification of the category
    if (angle<22.5):
        case = 3
    elif (angle<67.5):
        case = 0
    elif (angle<112.5):
        case = 1
    elif (angle<157.5):
        case = 2
    else:
        case = 3

    return case

def height_new(hgt):
    # id the categ
    if (hgt<5):
        hgt_cl = 1
    elif (hgt<10):
        hgt_cl = 2
    elif (hgt<15):
        hgt_cl = 3
    elif (hgt<20):
        hgt_cl = 4
    elif (hgt<25):
        hgt_cl = 5
    elif (hgt<30):
        hgt_cl = 6
    elif (hgt<35):
        hgt_cl = 7
    elif (hgt<40):
        hgt_cl = 8
    elif (hgt<50):
        hgt_cl = 9
    elif (hgt<60):
        hgt_cl = 10
    elif (hgt<80):
        hgt_cl = 11
    elif (hgt<100):
        hgt_cl = 12
    else:
        hgt_cl = 12

    return hgt_cl

def height_old(hgt):
    # id the categ
    if (hgt<5):
        hgt_cl = 0
    elif (hgt<10):
        hgt_cl = 1
    elif (hgt<15):
        hgt_cl = 2
    elif (hgt<20):
        hgt_cl = 3
    elif (hgt<25):
        hgt_cl = 4
    elif (hgt<30):
        hgt_cl = 5
    elif (hgt<35):
        hgt_cl = 6
    elif (hgt<40):
        hgt_cl = 7
    elif (hgt<50):
        hgt_cl = 8
    elif (hgt<60):
        hgt_cl = 9
    elif (hgt<80):
        hgt_cl = 10
    elif (hgt<100):
        hgt_cl = 11
    else:
        hgt_cl = 12

    return hgt_cl


