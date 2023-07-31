#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 19:00:59 2023

@author: konstantinos
"""

import numpy as np
import numba
from casters import THE_SMALL_CASTER
from time_extractor import days_since_distruption
from entropy_stream_identificator import entropy_id_vec
from astropy.coordinates import cartesian_to_spherical

#%% Specify what we are looking at

coord = 'r' # Can be r, theta, phi
quantity = 'vel' # Can be rho, vel
mass_weigh = True 
thing = coord + '-' + quantity
stream = False

Mbh = 4 # 6 or 4
if Mbh == 6:
    # Batch fixes so it's one call.
    part = 0
    fixes1 = np.arange(750, 820)
    fixes2 = np.arange(820, 900)
    fixes3 = np.arange(900, 970)
    fixes4 = np.arange(970,1008+1)
    fixes_list = [fixes1, fixes2, fixes3, fixes4]
    pre = 'tde_data/snap_'
if Mbh == 4:
    # Batch fixes so it's one call.
    part = 0
    fixes1 = np.arange(220, 320)
    fixes2 = np.arange(321, 412+1)
    fixes_list = [fixes1, fixes2]
    pre = 'tde_data2/snap_'
#%% Constants

G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2

# Need to convert velocity Rsol/time to km/s
Rsol_to_km = 6.957e5
converter_vel = Rsol_to_km / t 

# Need to convert density Msol/Rsol^3 to g/cm^3
Msol_to_g = 1.989e33
Rsol_to_cm = 6.957e10
converter_den = Msol_to_g / Rsol_to_cm**3

#%% Doing the thing

def THE_DECIDER(quantity, coord):
    if quantity == 'den':
        uncasted = np.load(pre + fix + '/Den_' + fix + '.npy')
        uncasted *= converter_den
        avg = False
        
    if quantity == 'vel':
        uncasted = V
        uncasted *= converter_vel
        avg = True
        
    if coord == 'r':
        old_grid = R
        r_start = 50
        r_stop = 80_000
        r_num = 100
        new_grid = np.linspace(r_start, r_stop, num = r_num)
        
    if coord == 'theta':
        old_grid = THETA
        theta_num = 100
        new_grid = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
            
    if coord == 'phi':
        old_grid = PHI
        phi_num = 100
        new_grid = np.linspace(-np.pi, np.pi, num = phi_num)
    
        
    return old_grid, new_grid, uncasted, avg


for fixes in fixes_list:
    part += 1
    if stream: 
        do_i_care_about_the_stream = [False, True]
    else:
        do_i_care_about_the_stream = [False]
    evolution = []
    evolution_stream = []
    for stream in do_i_care_about_the_stream:
        for fix in fixes:
            # Progress Check
            fix = str(fix)
            day = str(np.round(days_since_distruption(pre+fix+'/snap_'+fix+'.h5'),1))
            print('Day:', day)
            
            # CM Position Data
            X = np.load(pre + fix + '/CMx_' + fix + '.npy')
            Y = np.load(pre + fix + '/CMy_' + fix + '.npy')
            Z = np.load(pre + fix + '/CMz_' + fix + '.npy')
            
            # Convert to spherical
            R, THETA, PHI = cartesian_to_spherical(X,Y,Z)
            R = R.value
            THETA = THETA.value
            PHI = PHI.value
            PHI = (-1) * (PHI  - np.pi) # (0,2π) -> (-π,π)
            
            # Velocity Data
            Vx = np.load(pre + fix + '/Vx_' + fix + '.npy')
            Vy = np.load(pre + fix + '/Vy_' + fix + '.npy')
            Vz = np.load(pre + fix + '/Vz_' + fix + '.npy')
            V = np.sqrt(Vx**2 + Vy**2 + Vz**2)
                
            if stream:
            # Import Entropy, Mask the stream
                Entropy = np.load(pre + fix + '/Entropy_' + fix + '.npy')
                stream_mask = entropy_id_vec(Entropy) # Identify the stream
                stream_mask = (stream_mask - 1) * (-1) # 1 where there is NO stream
              
            	# We only care about unbound material
            Orbital = (0.5 * V**2 ) - 10**Mbh / (R-rg)
            unbound_mask = np.where(Orbital > 0, 1, 0) # returns 1 on unbound gas, 0 on bound
              
            # Specify new grid and decide the quantity to cast
            old_grid, new_grid, uncasted, avg = THE_DECIDER(quantity, coord)
            
            # Apply masks
            if stream:
                mask = np.multiply(stream_mask, unbound_mask)
                uncasted = np.multiply(uncasted, mask)
            else:
                uncasted = np.multiply(uncasted, unbound_mask)
              		
            # EVOKE THE CASTER
            Mass = np.load(pre + fix + '/Mass__' + fix + '.npy')
            casted = THE_SMALL_CASTER(new_grid, old_grid, uncasted, weights = Mass,
                                      avg=avg)
            casted = np.nan_to_num(casted)
              		          
            if stream:
                evolution_stream.append(casted)
            else:
                evolution.append(casted)
            
    np.save('tavg2/evo-M'+str(Mbh)+thing+'yes-'+str(part),evolution_stream)
    np.save('tavg2/evo-M'+str(Mbh)+thing+'no-'+str(part),evolution)
