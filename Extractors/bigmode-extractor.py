# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 21:08:41 2022

@author: Konstantinos
"""
#%% Imports

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.figsize'] = [10.0, 5.0]
plt.rcParams['figure.dpi'] = 300
from datetime import datetime
import h5py
import os

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2
#%% Get Energies

## File structure is
# box, cycle, time, mpi, rank0 ... rank99.
# This iterates over all the ranks
def calculator(f, key):
    Kin = []
    Pot = []
    Bound = []
    # for r
    x = f[key]['CMx']
    y = f[key]['CMy']
    z = f[key]['CMz']
    # for v
    vx = f[key]['Vx']
    vy = f[key]['Vy']
    vz = f[key]['Vz']

    for i in range(len(x)):
        v = np.linalg.norm([ vx[i], vy[i], vz[i]]) # * Rsol/t 
        KinE = 0.5*v**2 # [Rsol^2/t^2]
        # Reminder: G = 1, M_BH = 1e6 Msol
        r = np.linalg.norm([ x[i], y[i], z[i]]) # * Rsol
        # vesc = np.sqrt(2*Mbh/r)
        # Need Paczynski Witta Potential
        PotE =  -Mbh / (r-rg) # * G [Msol/Rsol]
        if np.abs(PotE) > KinE:
            Bound.append(1)
            Pot.append(PotE)
            Kin.append(KinE)
        else:
            Bound.append(0)
            Pot.append(PotE)
            Kin.append(KinE)
            
    return  Kin, Pot, Bound 

def combiner(fix, thing):
    
    # Load 2 parts
    part1 = np.load(fix + '/' + thing + '_' + fix + '_' + 'part1' + '.npy')
    part2 = np.load(fix + '/' + thing + '_' + fix + '_' + 'part2' + '.npy')

    # Concatenate            
    combined = np.concatenate((part1, part2))
    
    # Save
    np.save(fix + '/' + thing + '_' + fix + '_BIGMODE', combined)
    
    # Delete the constituent parts
    os.remove(fix + '/' + thing + '_' + fix + '_' + 'part1' + '.npy')
    os.remove(fix + '/' + thing + '_' + fix + '_' + 'part2' + '.npy')
    return
    
def extractor(fix, bigmode = False):
    '''
    Loads the file, extracts specific kinetic and potential energies 
    
    Parameters
    ----------
    fix : str, 
        hdf5 file suffix and prefix.
        
    bigmode : bool, default False
            Specifies whether to employ special procedures for large
            (>16 GB) files 
    Returns
    -------
    Kin : np.array, float64
        Specific Kinetic Energy 
    Pot : np.array, float64
        Specific Potential Energy
    Bound : np.array, bool
        True for bound material, false for unbound
        Bound means Pot > Kin
    
    '''
    # Timing start
    start_time = datetime.now()
    # Read File
    snapshot = fix+'/snap_'+fix+'.h5'
    f = h5py.File(snapshot, "r")
    # HDF5 are dicts, get the keys.
    keys = f.keys() 
    # List with keys that don't hold relevant data
    not_ranks = ['Box', 'Cycle', 'Time', 'mpi']
            
    # Bigmode for large datasets
    # all it does is that it splits the dataset in two, does the thing
    # for one half and saves output to a .npy
    # Then it does the same for the second half. Finally it combines the two 
    if bigmode:
        # Generate key list
        keys_list = []
        print(type(keys))
        for key in keys:
            if key not in not_ranks:
                keys_list.append(key)
        
        # Split key list in half
        halfpoint = len(keys_list)//2
        keys_1 = keys_list[:halfpoint]
        keys_2 = keys_list[halfpoint:]
        
        # First half lists
        Kin_1 = []
        Pot_1 = []
        Bound_1 = []
        
        # First half
        for key in keys_1:
            # Sanity Check & Timing
            printing_ranks = ['rank1','rank2','rank3', 'rank4','rank5',
                              'rank6','rank7','rank8','rank9']
            end_time = datetime.now()
            if key in printing_ranks:
                print(key)
                print('Duration: {}'.format(end_time - start_time))
            
            # Calculator
            Kin_temp, Pot_temp, Bound_temp = calculator(f, key)
            Kin_1 = Kin_1 + Kin_temp
            Pot_1 = Pot_1 + Pot_temp
            Bound_1 = Bound_1 + Bound_temp
    
        
        # Save to another file.
        np.save(fix+'/KinCM_'+fix+'_part1', Κin_1)   
        np.save(fix+'/PotCM_'+fix+'_part1', Pot_1) 
        np.save(fix+'/BoundCM_'+fix+'_part1', Bound_1)
        
        # Free some space.   
        del Kin_1, Pot_1, Bound_1 
        
        # Second half lists
        Kin_2 = []
        Pot_2 = []
        Bound_2 = []
        
        # First half
        for key in keys_2:
            # Sanity Check & Timing
            printing_ranks = ['rank1','rank2','rank3', 'rank4','rank5',
                              'rank6','rank7','rank8','rank9']
            end_time = datetime.now()
            if key in printing_ranks:
                print(key)
                print('Duration: {}'.format(end_time - start_time))
            
            # Call Calculator
            Kin_temp, Pot_temp, Bound_temp = calculator(f, key)
            Kin_2 = Kin_2 + Kin_temp
            Pot_2 = Pot_2 + Pot_temp
            Bound_2 = Bound_2 + Bound_temp
        
        # Close the file
        f.close()
        
        # Save to another file.
        np.save(fix+'/KinCM_'+fix+'_part2', Kin_2)   
        np.save(fix+'/PotCM_'+fix+'_part2', Pot_2) 
        np.save(fix+'/BoundCM_'+fix+'_part2', Bound_2)
        
        # Free some space.   
        del Kin_2, Pot_2, Bound_2
        
        # Combine the 2 parts
        combiner(fix, 'KinCM')
        combiner(fix, 'PotCM')
        combiner(fix, 'BoundCM')
        return
    #################################################
    # Small mode  
    # Use lists for clarity
    Kin = []
    Pot = []
    Bound = []
    for key in keys:
        if key in not_ranks:
            # Skip whatever is not a mpi rank
            continue
        else:
            # Sanity Check & Timing
            printing_ranks = ['rank1','rank2','rank4','rank7','rank9']
            end_time = datetime.now()
            if key in printing_ranks:
                print(key)
                print('Duration: {}'.format(end_time - start_time))
            # For some reason, having the collumns into variables is way faster.
            # for r
            x = f[key]['CMx']
            y = f[key]['CMy']
            z = f[key]['CMz']
            # for v
            vx = f[key]['Vx']
            vy = f[key]['Vy']
            vz = f[key]['Vz']
            # mass, not needed for specific energies
            # den = f[key]['Density']
            # vol = f[key]['Volume']
            # mass = den * vol
            for i in range(len(x)):
                v = np.linalg.norm([ vx[i], vy[i], vz[i]]) # * Rsol/t 
                KinE = 0.5*v**2 # [Rsol^2/t^2]
                # Reminder: G = 1, M_BH = 1e6 Msol
                r = np.linalg.norm([ x[i], y[i], z[i]]) # * Rsol
                # vesc = np.sqrt(2*Mbh/r)
                # Need Paczynski Witta Potential
                PotE =  -Mbh / (r-rg) # * G [Msol/Rsol]
                if np.abs(PotE) > KinE:
                    Bound.append(1)
                    Pot.append(PotE)
                    Kin.append(KinE)
                else:
                    Bound.append(0)
                    Pot.append(PotE)
                    Kin.append(KinE)
    # Close the file
    f.close()
    return


#%% Doing it
fix = '1008'
extractor(fix, bigmode=True)   
    
    
    
    
    
    
    
    
    
    
    
    
