#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 15:02:38 2023

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
    Mass = []
    
    # Iterate over ranks
    for key in keys:
        if key in not_ranks:
            # Skip whatever is not a mpi rank
            continue
        else:
            # For some reason, having the collumns into variables is way faster.
            vol_data = f[key]['Volume']
            den_data = f[key]['Density']
            for i in range(len(den_data)):
                Mass.append(vol_data[i] * den_data[i])

    # Close the file
    f.close()
    return Mass
#%%
# Change the current working directory
work_dir = input()
snap_num = input()
os.chdir(work_dir+'/')# NOTE: change backk to tde_data
fixes = np.arange(1,snap_num+1)
global_time = datetime.now()
for fix in fixes:
	# Timing start
	start_time = datetime.now()

	# Read file
	folder = 'snap_' + str(fix)+ '/'
	snapshot = folder + 'snap_' + str(fix) + '.h5'
	print(snapshot)
	Mass = extractor(snapshot)

	# Timing end
	end_time = datetime.now()
	print('Duration for this snapshot: {}'.format(end_time - start_time))
	print('Total time elapsed: {}'.format(end_time - global_time))   

	# Save to another file.
	np.save(folder + 'Mass__'+str(fix), Mass)   
