#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 19:31:03 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt

def gridder(start, stop, num, Rt):
    ''' Make a grid with more resolution near the BH'''
    num_fourths = int(num/4)    
    x1 = np.linspace(start, -3*Rt, num_fourths)
    x2 = np.linspace(-3*Rt, 3*Rt, 2*num_fourths)
    x3 = np.linspace(3*Rt, stop, num_fourths)
    xs = np.concatenate((x1, x2, x3))
    return xs

if __name__ == '__main__':
    x_start = -300
    x_stop = 50
    x_num = 520 # np.abs(x_start - x_stop)
    m = 4
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    xs = gridder(x_start, x_stop, x_num, Rt)