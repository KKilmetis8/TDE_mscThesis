#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 17:06:56 2023

@author: konstantinos
"""

import numpy as np
import h5py
from datetime import datetime
import os

#%% Get Densities

## File structure is
# box, cycle, time, mpi, rank0 ... rank99.
# This iterates over all the ranks



def extractor(filename):
    '''
    Loads the file, extracts X,Y,Z and Density. 
    
    Parameters
    ----------
    f : str, 
        hdf5 file name. Contains the data
    
    Returns
    -------
    X : np.array, float64
        X - coordinate
    Y : np.array, float64
        Y - coordinate.
    Z : np.array, float64
        Z - coordinate.
    Den : np.array, float64
        Density.
    
    '''
    # Timing start
    start_time = datetime.now()
    # Read File
    f = h5py.File(filename, "r")
    # HDF5 are dicts, get the keys.
    keys = f.keys() 
    # List to store the length of each rank
    lengths = []
    # List with keys that don't hold relevant data
    not_ranks = ['Box', 'Cycle', 'Time', 'mpi']
    
    for key in keys:
        if key in not_ranks:
            # Skip whatever is not a mpi rank
            continue
        else:
            # Store the length of the dataset
            lengths.append(len(f[key]['X']))
    
    # Use lists for clarity
    M = []
    
    # Iterate over ranks
    for key in keys:
        if key in not_ranks:
            # Skip whatever is not a mpi rank
            continue
        else:
            # Sanity Check
            print(key)
            # Timing
            end_time = datetime.now()
            print('Duration: {}'.format(end_time - start_time))
            # For some reason, having the collumns into variables is way faster.
            vol_data = f[key]['Volume']
            den_data = f[key]['Density']

            for i in range(len(den_data)):
                M.append(vol_data[i] * den_data[i])
    # Close the file
    f.close()
    return M
#%%

fixes = [217]
check = 's30'
sims = ['base-', check + '-']
for fix in fixes:
    for sim in sims:
        fix = str(fix)
        name = check + '/' + fix + '/' + sim + fix
        snapshot = name  + '.h5'
        print(snapshot)

        M = extractor(snapshot)  
        
        # Save to another file.
        np.save(name + 'Mass', M)   
            
