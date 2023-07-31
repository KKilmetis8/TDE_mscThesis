#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 17:43:02 2023

@author: konstantinos
"""

import numpy as np
import numba

def THE_CASTER_SPLIT(xs, X,
               ys, Y,
               Den, 
               weights=None, avg=False):
    '''
    Casts the density down to a smaller size vector

    Parameters
    ----------
    radii : arr,
        Array of radii we want to cast to.
    R : arr,
        Radius data.
    thetas: arr,
        Array of true anomalies we want to cast to.
    THETA: arr,
        True anomalies data
    Den: arr,
        Density data
        
    Returns
    -------
     density: arr
        Cast down version of density

    '''
    gridded_density = np.zeros((len(xs), len(ys)))
    gridded_weights = np.zeros((len(xs), len(ys)))
    counter = np.zeros((len(xs), len(ys)))
    current_progress = 0
    for i in range(len(X)):
        
        # Check how close the true R is to our radii
        diffs = np.abs(xs - X[i])
        
        # Get the index of the two closest cells in x
        diffs_sorted = np.argsort(diffs)
        idx_x1 = diffs_sorted[0]
        idx_x2 = diffs_sorted[1]
        
        # Same in y
        diffs = np.abs(ys - Y[i])
        diffs_sorted = np.argsort(diffs)
        idx_y1 = diffs_sorted[0]
        idx_y2 = diffs_sorted[1]
        
        # Get their Rs
        R = X[i]**2 + Y[i]**2
        rs = [ xs[idx_x1]**2 + ys[idx_y1]**2, 
               xs[idx_x1]**2 + ys[idx_y2]**2, 
               xs[idx_x2]**2 + ys[idx_y1]**2, 
               xs[idx_x2]**2 + ys[idx_y2]**2
               ]
        
        # Find how far away each is
        r_diffs = np.abs(R - rs)
        
        # Turn to a weight, giving precadence to the ones closest 
        inv_r_diffs = 1/r_diffs
        normalization_constant = np.sum(inv_r_diffs)
        splits = inv_r_diffs / normalization_constant
        # print(splits)
        
        # Add to grid hold for  now
        # counter[idx_x, idx_y] += 1
        
        # Add to grid | weighted mean
        if weights != None:
            
            gridded_density[idx_x1, idx_y1] += Den[i] * weights[i] * splits[0]
            gridded_weights[idx_x1, idx_y1] += weights[i] * splits[0]
            
            gridded_density[idx_x1, idx_y2] += Den[i] * weights[i] * splits[1]
            gridded_weights[idx_x1, idx_y2] += weights[i] * splits[1]
            
            gridded_density[idx_x2, idx_y1] += Den[i] * weights[i] * splits[2]
            gridded_weights[idx_x2, idx_y1] += weights[i] * splits[2]
            
            gridded_density[idx_x2, idx_y2] += Den[i] * weights[i] * splits[3]
            gridded_weights[idx_x2, idx_y2] += weights[i] * splits[3]
        else:
            gridded_density[idx_x1, idx_y1] += Den[i] * splits[0]
            gridded_density[idx_x1, idx_y2] += Den[i] * splits[1]
            gridded_density[idx_x2, idx_y1] += Den[i] * splits[2]
            gridded_density[idx_x2, idx_y2] += Den[i] * splits[3]

        
        # Progress check
        progress = int(np.round(i/len(X),1) * 100)
        if i % 100  == 0 and progress != current_progress:
            print('THE CASTER IS', progress, '% DONE')
            current_progress = progress
    # Normalize
    final_density = gridded_density
    if avg:
        final_density = np.divide(gridded_density,counter)
    if weights != None:
        final_density = np.divide(gridded_density, gridded_weights)
    return final_density

def THE_CASTER_SPLIT_2(xs, X,
               ys, Y,
               Den, 
               Vol, Z,
               weights=None, avg=False):
    '''
    Casts the density down to a smaller size vector

    Parameters
    ----------
    radii : arr,
        Array of radii we want to cast to.
    R : arr,
        Radius data.
    thetas: arr,
        Array of true anomalies we want to cast to.
    THETA: arr,
        True anomalies data
    Den: arr,
        Density data
        
    Returns
    -------
     density: arr
        Cast down version of density

    '''
    # Add only cells which |z| < cell_radius
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
    
    gridded_density = np.zeros((len(xs), len(ys)))
    gridded_weights = np.zeros((len(xs), len(ys)))
    counter = np.zeros((len(xs), len(ys)))
    current_progress = 0
    for i in range(len(X)):
        
        # Check how close the true R is to our radii
        diffs = np.abs(xs - X[i])
        
        # Get the index of the two closest cells in x
        diffs_sorted = np.argsort(diffs)
        idx_x1 = diffs_sorted[0]
        idx_x2 = diffs_sorted[1]
        
        # Same in y
        diffs = np.abs(ys - Y[i])
        diffs_sorted = np.argsort(diffs)
        idx_y1 = diffs_sorted[0]
        idx_y2 = diffs_sorted[1]
        
        # Get their Rs
        R = X[i]**2 + Y[i]**2
        rs = [ xs[idx_x1]**2 + ys[idx_y1]**2, 
               xs[idx_x1]**2 + ys[idx_y2]**2, 
               xs[idx_x2]**2 + ys[idx_y1]**2, 
               xs[idx_x2]**2 + ys[idx_y2]**2
               ]
        
        # Find how far away each is
        r_diffs = np.abs(R - rs)
        
        # Turn to a weight, giving precadence to the ones closest 
        inv_r_diffs = 1/r_diffs
        normalization_constant = np.sum(inv_r_diffs)
        splits = inv_r_diffs / normalization_constant
        # print(splits)
        
        # Add to grid hold for  now
        # counter[idx_x, idx_y] += 1
        
        # Add to grid | weighted mean
        if weights != None:
            
            gridded_density[idx_x1, idx_y1] += Den[i] * weights[i] * splits[0]
            gridded_weights[idx_x1, idx_y1] += weights[i] * splits[0]
            
            gridded_density[idx_x1, idx_y2] += Den[i] * weights[i] * splits[1]
            gridded_weights[idx_x1, idx_y2] += weights[i] * splits[1]
            
            gridded_density[idx_x2, idx_y1] += Den[i] * weights[i] * splits[2]
            gridded_weights[idx_x2, idx_y1] += weights[i] * splits[2]
            
            gridded_density[idx_x2, idx_y2] += Den[i] * weights[i] * splits[3]
            gridded_weights[idx_x2, idx_y2] += weights[i] * splits[3]
        else:
            gridded_density[idx_x1, idx_y1] += Den[i] * splits[0]
            gridded_density[idx_x1, idx_y2] += Den[i] * splits[1]
            gridded_density[idx_x2, idx_y1] += Den[i] * splits[2]
            gridded_density[idx_x2, idx_y2] += Den[i] * splits[3]

        
        # Progress check
        progress = int(np.round(i/len(X),1) * 100)
        if i % 100  == 0 and progress != current_progress:
            print('THE CASTER IS', progress, '% DONE')
            current_progress = progress
    # Normalize
    final_density = gridded_density
    if avg:
        final_density = np.divide(gridded_density,counter)
    if weights != None:
        final_density = np.divide(gridded_density, gridded_weights)
    return final_density