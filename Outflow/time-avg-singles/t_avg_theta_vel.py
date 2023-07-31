#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:55:45 2023

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:49:32 2023

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 19:05:16 2023

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 16:06:27 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8 , 4]
import numba
import matplotlib.colors as colors
import colorcet as cc # cooler colormaps
from calculators.casters import THE_SMALL_CASTER
from extractors.time_extractor import days_since_distruption
from astropy.coordinates import cartesian_to_spherical
from calculators.entropy_stream_identificator import entropy_id_vec
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
fixes = ['820', '881', '925']
theta_evolution = []
theta_evolution_stream = []
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
        
        R, THETA, PHI = cartesian_to_spherical(X,Y,Z)
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
        theta_num = 100
        thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
        
        # EVOKE
        v_cast = THE_SMALL_CASTER(thetas, THETA, Unbound_V, avg=True)
        v_cast = np.nan_to_num(v_cast)
        
        # Remove first and last element, they get smooshed
        v_cast = np.delete(v_cast, 0)
        v_cast = np.delete(v_cast, -1)
        
        if stream:
            theta_evolution_stream.append(v_cast)
        else:
            theta_evolution.append(v_cast)
#%% Time average
@numba.njit
def time_averager(evolution):
    time_averaged = np.zeros(len(evolution[0]))
    for i in range(len(evolution[0])): # Loop over radii
        for j in range(len(evolution)): # Loop over time
            time_averaged[i] += evolution[j][i]
        time_averaged[i] /= len(evolution)
    return time_averaged

time_averaged_theta = time_averager(theta_evolution)
time_averaged_theta_stream = time_averager(theta_evolution_stream)

# Plotting
thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
thetas = np.delete(thetas, 0)
thetas = np.delete(thetas, -1)
fig, ax = plt.subplots()
ax.plot(thetas, time_averaged_theta, 
         color = 'navy', label = 'Yes Stream')
ax.plot(thetas, time_averaged_theta_stream,
        color = 'darkorange', linestyle='dashed', label = 'No Stream')
ax.set_title(r'Time averaged $\textbf{outflow}$ velocity', fontsize = 18)
ax.set_xlabel(r'$\theta$ [rad] (Azimuthal Angle)', fontsize = 14)
ax.set_ylabel(r'$ |V| \left[ \frac{\mathrm{km}}{\mathrm{s}} \right] $',fontsize = 14)
ax.set_yscale('log')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.7, 0.8, 'This uses only 3 snapshots!, \n It will be deployed at scale.',
         transform=ax.transAxes, bbox = props)
plt.legend(loc='upper left')
plt.grid()
    
    