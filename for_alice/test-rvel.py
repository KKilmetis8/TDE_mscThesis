#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 19:14:11 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Serif'
plt.rcParams['figure.figsize'] = [8 , 4]
import numba
import matplotlib.colors as colors
from casters import THE_SMALL_CASTER
from time_extractor import days_since_distruption
from entropy_stream_identificator import entropy_id_vec
#%% Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2

# Need to convert density Msol/Rsol^2 to g/cm
Rsol_to_km = 6.957e5
converter = Rsol_to_km / t

#%% Doing the thing
do_i_care_about_the_stream = [False, True]
fixes = ['820'] #, '881', '925']
rho_evolution = []
rho_evolution_stream = []
pre = 'tde_data/snap_'
for stream in do_i_care_about_the_stream:
    for fix in fixes:
        # Progress Check
        day = str(np.round(days_since_distruption(pre+ fix+'/snap_'+fix+'.h5'),1))
        print('Day:', day)
        # CM Position Data
        X = np.load(pre + fix + '/CMx_' + fix + '.npy')
        Y = np.load(pre + fix + '/CMy_' + fix + '.npy')
        Z = np.load(pre + fix + '/CMz_' + fix + '.npy')
        R = np.sqrt(X**2 + Y**2 + Z**2)
        del X, Y, Z
        # Velocity Data
        Vx = np.load(pre + fix +'/Vx_' + fix + '.npy')
        Vy = np.load(pre + fix +'/Vy_' + fix + '.npy')
        Vz = np.load(pre + fix +'/Vz_' + fix + '.npy')
        V = np.sqrt(Vx**2 + Vy**2 + Vz**2)
        Mass = np.load(pre + fix + '/Mass__' + fix + '.npy')
        del Vx, Vy, Vz
        
        if stream:
            # Import Entropy, Mask the stream
            Entropy = np.load(pre + fix + '/Entropy_' + fix + '.npy')
            stream_mask = entropy_id_vec(Entropy) # Identify the stream
            stream_mask = (stream_mask - 1) * (-1) # 1 where there is NO stream
        
        # We only care about unbound material
        Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
        unbound_mask = np.where(Orbital > 0, 1, 0) # returns 1 on unbound gas, 0 on bound
        V *= converter # Convert to km/s
        
        # Combine the masks
        if stream:
            mask = np.multiply(stream_mask, unbound_mask)
            Unbound_V = np.multiply(V, mask)
        else:
            Unbound_V = np.multiply(V, unbound_mask)
            
        # Specify new grid:
        r_start = 50
        r_stop = 80_000
        r_num = 100 # np.abs(x_start - x_stop)
        radii = np.linspace(r_start, r_stop, num = r_num)
        
        # EVOKE
        v_cast = THE_SMALL_CASTER(radii, R, Unbound_V, 
                                  weights = Mass, avg = True) # cant add V
        v_cast = np.nan_to_num(v_cast)
        if stream:
            rho_evolution_stream.append(v_cast)
        else:
            rho_evolution.append(v_cast)

np.save('test',rho_evolution)
np.save('test_s',rho_evolution_stream)
#%% Time average
# @numba.njit
# def time_averager(evolution):
#     time_averaged = np.zeros(len(evolution[0]))
#     for i in range(len(evolution[0])): # Loop over radii
#         for j in range(len(evolution)): # Loop over time
#             time_averaged[i] += evolution[j][i]
#         time_averaged[i] /= len(evolution)
#     return time_averaged

# time_averaged_rho = time_averager(rho_evolution)
# time_averaged_rho_stream = time_averager(rho_evolution_stream)

# Plotting
rho_evolution = np.load('test2-2.npy')
rho_evolution_stream = np.load('test2_s-2.npy')
r_start = 50
r_stop = 80_000
r_num = 100 # np.abs(x_start - x_stop)
radii = np.linspace(r_start, r_stop, num = r_num)
# radii = np.delete(radii, 0)
# radii = np.delete(radii, -1)

fig, ax = plt.subplots()
ax.plot(radii, rho_evolution.T, 
         color = 'purple', label='Yes Stream')
ax.plot(radii, rho_evolution_stream.T, 
         color = 'goldenrod', linestyle='--', label='No Stream')
ax.set_title(r'Time averaged, Mass-Weighed $\textbf{outflow}$ velocity', fontsize = 18)
ax.set_xlabel(r'Radial Distance $\left[ R_{\odot} \right]$ ', fontsize = 14)
ax.set_ylabel(r'$ |V| \left[ \frac{\mathrm{km}}{\mathrm{s}} \right] $',fontsize = 14)
ax.set_yscale('log')
ax.set_ylim(1e-2,1e6)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.4, 0.2, 'This is only snapshot 36.6! \n It will be deployed at scale.',
         transform=ax.transAxes, bbox = props)
plt.grid()
plt.legend()
plt.savefig('test')
    