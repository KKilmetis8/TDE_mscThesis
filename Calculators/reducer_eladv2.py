#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 4 21:14:07 2023

@author: konstantinos
"""

import numpy as np
import numba

def THE_REDUCER_2(xs, X,
                ys, Y,
                Den,
                Vol, Z):
    
    # Add only cells which |z| < 2*cell_radius
    cell_radius = np.power(0.25 * 3 * Vol /np.pi, 1/3)
    new_X = []
    new_Y = []
    for i in range(len(Z)):
        diff = np.abs(Z[i]) - 2*cell_radius[i]
        if diff < 0:
            new_X.append(X[i])
            new_Y.append(Y[i])
    X = new_X
    Y = new_Y
    
    # Make pixels
    pixels = []
    for x in xs:
        for y in ys:
            pixels.append( (x,y) )

    ungridded_den = np.zeros( (len(xs) * len(ys)) )
    for i in range(len(pixels)):
        
        # Find closest
        diffs = (pixels[i][0] - X)**2 + (pixels[i][1] - Y)**2 
        idx = np.argmin(diffs)
        
        # Add
        ungridded_den[i] = Den[idx]
        
        # Progress
        progress = 100 * i/len(pixels)
        if progress - int(progress) == 0:
            print(int(progress), ' % done')

    # Reshape
    gridded_density = np.reshape(ungridded_den, (len(xs), len(ys)))
            
    return gridded_density
