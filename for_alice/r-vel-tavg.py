#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 13:37:05 2023

@author: konstantinos
"""

import numpy as np
import numba
import matplotlib.colors as colors
from casters import THE_SMALL_CASTER
from time_extractor import days_since_distruption
from entropy_stream_identificator import entropy_id_vec
#%% Constants
# Need to convert density Msol/Rsol^2 to g/cm
Msol_to_g = 1.989e33
Rsol_to_cm = 6.957e10
converter = Msol_to_g / Rsol_to_cm**2

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2

#%% Doing the thing
pre = 'r-vel-'
do_i_care_about_the_stream = [False, True]
# fixes = np.arange(820,880+1) # 36-45
# name = '36-45'
# fixes = np.arange(881, 953+1) # 45-55
# name = '45-55'
fixes = np.arange(954,1008+1) # 55-end
name = '55+'
rho_evolution = []
rho_evolution_stream = []
for stream in do_i_care_about_the_stream: 
    for fix in fixes:
  	# Progress Check
        fix = str(fix)
        day = str(np.round(days_since_distruption('tde_data/snap_' +fix+'/snap_'+fix+'.h5'),1))
        print('Day:', day)
        # CM Position Data
        X = np.load('tde_data/snap_' + fix + '/CMx_' + fix + '.npy')
        Y = np.load('tde_data/snap_' + fix + '/CMy_' + fix + '.npy')
        Z = np.load('tde_data/snap_' + fix + '/CMz_' + fix + '.npy')
        R = np.sqrt(X**2 + Y**2 + Z**2)
        # Calc true anomalies
        del X, Y, Z
        # Velocity Data
        Vx = np.load('tde_data/snap_'+ fix + '/Vx_' + fix + '.npy')
        Vy = np.load('tde_data/snap_'+ fix + '/Vy_' + fix + '.npy')
        Vz = np.load('tde_data/snap_'+ fix + '/Vz_' + fix + '.npy')
        V = np.sqrt(Vx**2 + Vy**2 + Vz**2)
        del Vx, Vy, Vz
          
        if stream:
        # Import Entropy, Mask the stream
            Entropy = np.load('tde_data/snap_' + fix + '/Entropy_' + fix + '.npy')
            stream_mask = entropy_id_vec(Entropy) # Identify the stream
            stream_mask = (stream_mask - 1) * (-1) # 1 where there is NO stream
          
          		# We only care about unbound material
        Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
        unbound_mask = np.where(Orbital > 0, 1, 0) # returns 1 on unbound gas, 0 on bound
        V *= converter # Convert to km/s 
          		# Combine the masks
        if stream:
            mask = np.multiply(stream_mask, unbound_mask)
            Unbound_v = np.multiply(V, mask)
        else:
            Unbound_v = np.multiply(V, unbound_mask)
          
        # Specify new grid:
        r_start = 50
        r_stop = 80_000
        r_num = 100 # np.abs(x_start - x_stop)
        radii = np.linspace(r_start, r_stop, num = r_num)
          		
        # EVOKE
        v_cast = THE_SMALL_CASTER(radii, R, Unbound_v)
        v_cast = np.nan_to_num(v_cast)
          		
        # Remove first and last element, they get smooshed
        v_cast = np.delete(v_cast, 0)
        v_cast = np.delete(v_cast, -1)
          
        if stream:
            rho_evolution_stream.append(v_cast)
        else:
            rho_evolution.append(v_cast)
#%% Time average
@numba.njit
def time_averager(evolution):
    time_averaged = np.zeros(len(evolution[0]))
    for i in range(len(evolution[0])): # Loop over radii
        for j in range(len(evolution)): # Loop over time
            time_averaged[i] += evolution[j][i]
        time_averaged[i] /= len(evolution)
    return time_averaged

time_averaged_rho = time_averager(rho_evolution)
time_averaged_rho_stream = time_averager(rho_evolution_stream)

np.save('tavg/'+pre+'yes'+name,time_averaged_rho)
np.save('tavg/'+pre+'no'+name,time_averaged_rho_stream)