#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 19:05:17 2023

@author: konstantinos
"""

import numpy as np
import numba

@numba.njit
def averager(xs):
    '''xs is a tuple containing the same physical quantity from
       a few different snapshots
    '''
    
    # Holds the result
    holder = np.zeros_like(xs[0]) 
    
    # Add all the snapshots
    for x in xs:
        holder = np.sum(holder, x)
        
    # Divide with the number of snapshots
    holder = np.multiply(holder, 1/len(xs))
    return holder