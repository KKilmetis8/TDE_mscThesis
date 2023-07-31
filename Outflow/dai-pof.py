#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 11:14:01 2023


@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# Pretty plots 
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8 , 4]
plt.rcParams['axes.facecolor']='whitesmoke'
import numba
# Custom
from calculators.casters import THE_CASTER
from extractors.time_extractor import days_since_distruption
from astropy.coordinates import cartesian_to_spherical
from calculators.entropy_stream_identificator import entropy_id_vec

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2

# Need to convert velocity Rsol/time to km/s
Rsol_to_km = 6.957e5
converter = Rsol_to_km / t

#%%
    
do_i_care_about_the_stream = [False, True]
fixes = ['820', '881']
dai_3 = []
dai_3_stream = []
for stream in do_i_care_about_the_stream:
    for fix in fixes:
        # Progress Check
        day = str(np.round(days_since_distruption(fix+'/snap_'+fix+'.h5'),1))
        print('Day:', day)
        # CM Position Data
        X = np.load(fix + '/CMx_' + fix + '.npy')
        Y = np.load(fix + '/CMy_' + fix + '.npy')
        Z = np.load(fix + '/CMz_' + fix + '.npy')
        R = np.sqrt(X**2 + Y**2 + Z**2)
        # Velocity Data
        Vx = np.load(fix + '/Vx_' + fix + '.npy')
        Vy = np.load(fix + '/Vy_' + fix + '.npy')
        Vz = np.load(fix + '/Vz_' + fix + '.npy')
        V = np.sqrt(Vx**2 + Vy**2 + Vz**2)
        del Vx, Vy, Vz
        Mass = np.load(fix +'/Mass_' + fix + '.npy')
        R, THETA, PHI = cartesian_to_spherical(X,Y,Z)
        R = R.value
        THETA = THETA.value
        
        if stream:
            # Import Entropy, Mask the stream
            Entropy = np.load(fix + '/Entropy_' + fix + '.npy')
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
        theta_num = 7
        thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
        
        # EVOKE
        v_cast = THE_CASTER(radii, R, thetas, THETA, Unbound_V, 
                            weights = Mass, avg=True)
        v_cast = np.nan_to_num(v_cast)
        
        # Remove first and last element, they get smooshed
        # v_cast = np.delete(v_cast, 0)
        # v_cast = np.delete(v_cast, -1)
        
        if stream:
            dai_3_stream.append(v_cast)
        else:
            dai_3.append(v_cast)
#%% Plotting
@numba.njit
def time_averager(evolution):
    time_averaged = np.zeros( np.shape(evolution[0]) )
    for i in range(len(evolution)): # Loop over time
        time_averaged = np.add(time_averaged, evolution[i])
    time_averaged = np.multiply(time_averaged, 1/len(evolution))
    return time_averaged

tavg = time_averager(dai_3)
tavg_stream = time_averager(dai_3_stream)

#%% Plotting
r_start = 50
r_stop = 80_000
r_num = 100 # np.abs(x_start - x_stop)
radii = np.linspace(r_start, r_stop, num = r_num)
theta_num = 7
thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
colarr = ['red', 'darkorange', 'gold',  'lightgreen', 'green', 'blue', 'violet']
line = 'solid'
fig, ax = plt.subplots(1,2, tight_layout=True)
for i in range(theta_num):
#     if i == 3:
#         line = 'dashed'
#         continue
    label = thetas[i] * 180 / np.pi
    label = str(np.round(label,0))
    ax[0].plot(radii, tavg[:,i], color = colarr[i], linestyle = line,
               label = label)
    ax[1].plot(radii, tavg_stream[:,i], color=colarr[i], linestyle = line,)
ax[0].grid()
ax[0].set_title('Stream Included')
ax[1].grid()
ax[1].set_title('Stream Excluded')

ax[0].set_xlabel(r'Radial Coordinate [$R_{\odot}$]')
ax[0].set_ylabel('Velocity [km/s]')
ax[1].set_xlabel(r'Radial Coordinate [$R_{\odot}$]')
#ax[1].set_ylabel('Velocity [km/s]')

ax[0].set_yscale('log')
ax[1].set_yscale('log')
# ax[0].set_ylim(2,2.5e4)
# ax[1].set_ylim(2,2.5e4)
fig.suptitle(r'Binned look at \textbf{outflow} velocity', fontsize = 16)
fig.legend(loc=7, title = r'$\theta$ bins')
fig.tight_layout()
fig.subplots_adjust(right=0.88)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
fig.text(0.8, 1.05, 'This uses only 2 snapshots! \n It will be deployed at scale.',
         transform=ax[0].transAxes, bbox = props)