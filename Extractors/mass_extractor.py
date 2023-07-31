# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 15:35:34 2022

@author: Konstantinos
"""

import numpy as np
import h5py
from datetime import datetime

#%% Get Mass

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
    M : np.array, float64
        Mass
    
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
    
    # Get the total length of each array.
    # total_length = np.sum(lengths)
    
    # Use numpy for speed
    # Den = np.zeros(total_length) 
    # X = np.zeros(total_length)
    # Y = np.zeros(total_length)
    # Z = np.zeros(total_length)
    # curr_len = 0
    
    # Use lists for clarity
    Mass = []
    Vol = []
    
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
            v_data = f[key]['Volume']
            den_data = f[key]['Density']
            for i in range(len(den_data)):
                Vol.append(v_data[i])
                Mass.append(v_data[i] * den_data[i])

    # Close the file
    f.close()
    return Mass, Vol
    
fixes = [445]
check = 'hr2'
# sims = [check + '-'] 
sims = ['base-', check + '-']
for fix in fixes:
    for sim in sims:
        fix = str(fix)
        name = 'convergence/' + fix + '/' + sim + fix
        snapshot = name  + '.h5'
        print(snapshot)

        Mass, Vol = extractor(snapshot)  
        
        # Save to another file.
        np.save(name + 'Mass', Mass)
        np.save(name + 'Volume', Vol)
            


    
    
    
    
    
    
    
    
    
    
    
    
    
    