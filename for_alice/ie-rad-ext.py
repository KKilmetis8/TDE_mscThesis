#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 15:04:03 2023

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
    T = []
    Y = []
    Z = []
    Den = []
    
    # Iterate over ranks
    for key in keys:
        if key in not_ranks:
            # Skip whatever is not a mpi rank
            continue
        else:
            # For some reason, having the collumns into variables is way faster.
            T_data = f[key]['InternalEnergy']
            Z_data = f[key]['tracers']['ZRadEnergy']
            for i in range(len(T_data)):
                T.append(T_data[i])
                Z.append(Z_data[i])

    # Close the file
    f.close()
    return T, Z
#%%
# Change the current working directory
os.chdir('tde_data/')
fixes = np.arange(750,1008+1)
global_time = datetime.now()
for fix in fixes:
	# Timing start
	start_time = datetime.now()

	# Read file
	folder = 'snap_' + str(fix)+ '/'
	snapshot = folder + 'snap_' + str(fix) + '.h5'
	print(snapshot)
	T, rad = extractor(snapshot)
	np.save(folder + 'IE_'+str(fix), T)   
	np.save(folder + 'Rad_'+str(fix), rad)   

	# Timing end
	end_time = datetime.now()
	print('Duration for this snapshot: {}'.format(end_time - start_time))
	print('Total time elapsed: {}'.format(end_time - global_time))   
