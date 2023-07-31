#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 21:14:07 2023

@author: konstantinos
"""

import numpy as np
import numba
import scipy.stats as stats
import time

def THE_REDUCER_3(xs, X,
                  ys, Y, Den):
    grid = [[], []]
    
    # Make 1-D grid
    for i, x in enumerate(xs):
        for j, y in enumerate(ys):
            grid[0].append(x)
            grid[1].append(y)

    # KDE
    values = np.vstack([X, Y])
    kernel = stats.gaussian_kde(values, weights = Den)
    not_gridded_density = kernel(grid)
    
    # Turn to 2-D
    gridded_density = np.zeros((len(xs), len(ys)))
    k = 0
    for i in range(len(xs)):
        for j in range(len(ys)):
            gridded_density[i,j] = not_gridded_density[k]
            k += 1
    
    return gridded_density

def THE_FILTERER(X, Y, Z, Den, Vol):
    # Add only cells which |z| < 2*cell_radius
    if type(Vol) != type(None):
        cell_radius = np.power(0.25 * 3 * Vol /np.pi, 1/3)
        new_X = []
        new_Y = []
        new_Z = []
        new_Den = []
        
        for i in range(len(Z)):
            diff = np.abs(Z[i]) - 2*cell_radius[i]
            if diff < 0:
                new_X.append(X[i])
                new_Y.append(Y[i])
                new_Z.append(Z[i])
                new_Den.append(Den[i])
        return new_X, new_Y, new_Z, new_Den
    
def THE_PIXELATOR(xs, ys, zs):
    # Make pixels
    pixels = []
    for x in xs:
        for y in ys:
            for z in zs:
                pixels.append( (x,y,z) )
    return pixels

@numba.njit
def THE_SUBTRACTOR(x, X, y, Y, z, Z):
    x_diff = np.power( np.add(x,X), 2)
    y_diff = np.power( np.add(y,Y), 2)
    z_diff = np.power( np.add(z,Z), 2)
    diffs = np.add(x_diff, y_diff) 
    diffs = np.add(diffs, z_diff)

    return diffs

def THE_REDUCER_2(pixels, shape, X, Y, Z, Den):
    # Make them negative
    X = np.array(X)
    X = np.multiply(X, -1)
    Y = np.array(Y)
    Y = np.multiply(Y, -1)
    Z = np.array(Z)
    Z = np.multiply(Z, -1)
    
    # Initialize lengths and holders
    old_progress = -1
    total_pixels = len(pixels)
    ungridded_den = np.zeros((total_pixels)) 
    for i in range(total_pixels):
        
        # Find closest 75 is probelem line for numbas and the time bottleneck
        # t = time.time()
        x = pixels[i][0]
        y = pixels[i][1]
        z = pixels[i][2]
        diffs = THE_SUBTRACTOR(x, X, y, Y, z, Z)
        idx = np.argmin(diffs)
        # t1 = time.time()

        # t2 = time.time()
        # print('Diffs: ',t1 - t)
        # print('Argmin: ', t2 - t1)
        
        # Add
        ungridded_den[i] = Den[idx]
        
        # Progress
        progress = int(100 * i/total_pixels)
        if progress % 10 == 0 and progress != old_progress:
            old_progress = progress
            print(int(progress), '% done')
    # Reshape
    gridded_density = np.reshape(ungridded_den, shape)            
    return gridded_density #gridded_density

def THE_PROJECTOR(den_cast, third_ax, ax):
    shape  = np.shape(den_cast)
    if ax == 'XY':
        den_project = np.zeros( (shape[0], shape[1]) )
        for i in range(len(third_ax)):
            den_project = np.add(den_project, den_cast[...,i])
    
        return den_project
    if ax == 'XZ':
        den_project = np.zeros( (shape[0], shape[2]) )
        for i in range(len(third_ax)):
            den_project = np.add(den_project, den_cast[:,i,:])
    
        return den_project
    if ax == 'YZ':
        den_project = np.zeros( (shape[1], shape[2]) )
        for i in range(len(third_ax)):
            den_project = np.add(den_project, den_cast[i,...])
    
        return den_project
