# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 18:32:48 2022

@author: Konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.figsize'] = [10.0, 5.0]
plt.rcParams['figure.dpi'] = 300
import pandas as pd
import colorcet as cc
from datetime import datetime

#%% Data Import

fix = '820' # It's both a prefix and a suffix so they cancel out!
X = np.load(fix + '/X_' + fix + '.npy')
Y = np.load(fix + '/Y_' + fix + '.npy')
Z = np.load(fix + '/Z_' + fix + '.npy')
Mass = np.load(fix + '/Mass_' + fix + '.npy')

#%% Aux. Functions

# Function to round to closest grid point

from bisect import bisect_left

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    
    NOTE: Taken from Stackoverflow wholesale
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before
    
def fix_ticks(a_range, b_range):
    '''
    EXTREMELY sloppy. Converts ticks from index number to solar radii
    requires a plot already to be present
    
    HAS to be called from within projector. This exists to make the 
    projector function more readable

    Parameters
    ----------
    a_range : array-like,
        Cell coordinates of a-axis.
    b_range : array-like,
        Cell coordinates of b-axis.

    Returns
    -------
    It fixes the ticks.

    '''
    # # Convert ticks to solar radii
    # # Get current ticks
    xlocs, xlabels = plt.xticks() 
    # turn them to ints to use as indices
    xlocs = xlocs.astype(np.int16)
    # get new labels
    new_x_labels = a_range[xlocs[1:-1]]
    # turn to int for aesthetics
    new_x_labels = new_x_labels.astype(np.int16)
    # give ticks
    plt.xticks(xlocs[1:-1]
                ,labels = new_x_labels)
               
    # Get current ticks
    ylocs, ylabels = plt.yticks() 
    # turn them to ints to use as indices
    ylocs = ylocs.astype(np.int16)
    # get new labels
    new_y_labels = b_range[ylocs[1:-1]]
    # turn to int for aesthetics
    new_y_labels = new_y_labels.astype(np.int16)
    # give ticks
    plt.yticks(ylocs[1:-1]
                ,labels = new_y_labels)

#%% Create Uniform Grid

def projector(axis_a, axis_b, Mass, max_a, min_a, max_b, min_b,
              point_num_a = 0,
              point_num_b = 0):
    '''
    Plots a projection to two axes by:
        1. Creating a uniform grid
        2. Adding mass to each cell
        3. Dividing by cell length
    NOTE: Can definetly be optimized by
        1. arrays instead of dataframes
        2. numbas
        3. sort by one axis, after going out of its maximum bound break the
           for loop.
           
    Parameters
    ----------
    axis_a : 1-D array-like, 
        first of the axes to project to.
    axis_b : 1-D array-like, 
        second of the axes to project to.
    Mass : 1-D array-like,
        Mass.
    max_a : int,
        maximum value of axis a to plot
    min_a : int,
        minimum value of axis a to plot
    max_b : int,
        maximum value of axis b to plot
    min_b : int,
        minimum value of axis b to plot.
    point_num_a : int, optional
        How many pixels the image should have in the horizontal axis.
        The default is |max_a - min_a|.
    point_num_b : int, optional
        How many pixels the image should have in the vertical
        axis.
        The default is |max_b - min_b|.
        

    Returns
    -------
    array : 2-D array
        Contains the projected density for the given axes.

    '''
    # Generates the grid length for each axis if nothing was inputted
    if point_num_a == 0 :
        point_num_a = np.abs(max_a - min_a) 
    if  point_num_b == 0:
        point_num_b = np.abs(max_b - min_b)
        
    # DF to np, only if I input a df series.
    if isinstance(axis_a, pd.core.series.Series):
        axis_a = axis_a.to_numpy()
        axis_b = axis_b.to_numpy()
        Mass = Mass.to_numpy()
        
    # Need to convert Msol/Rsol^2 to g/cm
    Msol_to_g = 1.989e33
    Rsol_to_cm = 6.957e10
    converter = Msol_to_g / Rsol_to_cm**2
    
    # Make ranges
    a_range = np.linspace(min_a, max_a, num = point_num_a)
    a_len = len(a_range) / (a_range[-1] - a_range[0])
    b_range = np.linspace(min_b, max_b, num = point_num_b)
    b_len = len(b_range) / (b_range[-1] - b_range[0])
    
    # Make grid
    grid = np.zeros( (point_num_a, point_num_b)) # 2D n x n array to store mass
    # Store everything in a DF
    # Flip needed to start from minimum
    df = pd.DataFrame(data = grid, index = a_range, columns = np.flip(b_range))
    # Assign mass to grid points
    # can definently be done faster if I used only np
    old_done = 0
    for i in range(len(Mass)):
        # Ensures relevant axis-a range
        if axis_a[i] < min_a or axis_a[i] > max_a:
            continue
        # Ensures relevant axis-b range
        elif axis_b[i] < min_b or axis_b[i] > max_b:
            continue
        else:
            grid_a = take_closest(a_range, axis_a[i] / a_len)
            grid_b = take_closest(b_range, axis_b[i] / b_len)
            
            df.loc[grid_a][grid_b] += Mass[i]
        
        # Progress Check
        done = int(100*i/len(Mass))
        if done>old_done:
            old_done = done    
            print( done, ' percent through the file')
    
    # Wrangle the data a bit
    
    # Convert to array for quicker plotting
    array = df.to_numpy()
    # Convert to log density by /dy dz
    array = np.log10(array * converter / (a_len*b_len) )
    
    array = np.nan_to_num(array.T, copy=True, nan=0.0, posinf=0.0,
                          neginf=0)
    # Renormalize Colors
    array[array > 5] = 5
    array[array < 0.3] = 0

    # Plot
    # Colormaps
    fire = cc.cm.fire
    rainbow = cc.cm.rainbow
    
    fig, ax = plt.subplots()
    img = ax.imshow(array, cmap=fire)
    plt.colorbar(img)
    
    # Fix ticks
    fix_ticks(a_range, b_range)
    return ax

#%% Mass
# Let's input only the first N data points, to help with time
simdata = pd.DataFrame({'X':X,
                   'Y':Y,
                   'Z':Z,
                   'Mass':Mass})
# Sort by mass
simdata.sort_values(by='Mass', ascending=False, inplace=True)
# Project for the N most massive points
start_time = datetime.now()
end_point = len(Mass)
xyproj = projector(simdata['X'][:end_point], 
                   simdata['Y'][:end_point], 
                   simdata['Mass'][:end_point], 
                  min_a = -15_000, max_a = 2000,
                  min_b = -4_000, max_b = 4_000,
                  )
end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))
xyproj.set_title(r'Logarithmic Density -$\log_{10}(\rho)$ [g/cm$^2$]- Projection for Snapshot 881 in the XY Plane')
xyproj.set_xlabel('X - Coordinate [$R_{\odot}$]')
xyproj.set_ylabel('Y - Coordinate [$R_{\odot}$]')
plt.show(xyproj)
