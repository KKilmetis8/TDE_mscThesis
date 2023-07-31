#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 18:26:57 2023

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 18:02:42 2023

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
plt.rcParams['figure.figsize'] = [8 , 4]
plt.rcParams['font.family'] = 'Times New Roman'
import numba
import matplotlib.colors as colors
import colorcet as cc # cooler colormaps
from calculators.casters import THE_SMALL_CASTER
from extractors.time_extractor import days_since_distruption
from astropy.coordinates import cartesian_to_spherical
from calculators.entropy_stream_identificator import entropy_id_vec
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
do_i_care_about_the_stream = [False, True]
fixes = ['820', '881', '925']
phi_evolution = []
phi_evolution_stream = []
for stream in do_i_care_about_the_stream:
    for fix in fixes:
        # Progress Check
        day = str(np.round(days_since_distruption(fix+'/snap_'+fix+'.h5'),1))
        print('Day:', day)
        # CM Position Data
        X = np.load(fix + '/CMx_' + fix + '.npy')
        Y = np.load(fix + '/CMy_' + fix + '.npy')
        Z = np.load(fix + '/CMz_' + fix + '.npy')
        
        R, THETA, PHI = cartesian_to_spherical(X,Y,Z)
        # Store it in an array of R, θ, φ, log(ρ)
        R = R.value 
        # THETA = THETA.value
        PHI = (-1) * (PHI.value  - np.pi) # (0,2π) -> (-π,π)
        del X, Y, Z
        # Velocity Data
        Vx = np.load(fix + '/Vx_' + fix + '.npy')
        Vy = np.load(fix + '/Vy_' + fix + '.npy')
        Vz = np.load(fix + '/Vz_' + fix + '.npy')
        V = np.sqrt(Vx**2 + Vy**2 + Vz**2)
        del Vx, Vy, Vz
        
        # Import Density
        Den = np.load(fix + '/Den_' + fix + '.npy')
        Den *=  converter
                
        if stream:
            # Import Entropy, Mask the stream
            Entropy = np.load(fix + '/Entropy_' + fix + '.npy')
            stream_mask = entropy_id_vec(Entropy) # Identify the stream
            stream_mask = (stream_mask - 1) * (-1) # 1 where there is NO stream
            
        # We only care about unbound material
        Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
        unbound_mask = np.where(Orbital > 0, 1, 0) # returns 1 on unbound gas, 0 on bound
        
        # Combine the masks
        if stream:
            mask = np.multiply(stream_mask, unbound_mask)
            Unbound_Den = np.multiply(Den, mask)
        else:
            Unbound_Den = np.multiply(Den, unbound_mask)
        
        # Specify new grid:
        phi_num = 100
        phis = np.linspace(-np.pi, np.pi, num = phi_num)
        
        # EVOKE
        den_cast = THE_SMALL_CASTER(phis, PHI, Unbound_Den)
        den_cast = np.nan_to_num(den_cast)
        
        # Remove first and last element, they get smooshed
        den_cast = np.delete(den_cast, 0)
        den_cast = np.delete(den_cast, -1)
        if stream:
            phi_evolution_stream.append(den_cast)
        else:
            phi_evolution.append(den_cast)
#%% Time average
@numba.njit
def time_averager(evolution):
    time_averaged = np.zeros(len(evolution[0]))
    for i in range(len(evolution[0])): # Loop over radii
        for j in range(len(evolution)): # Loop over time
            time_averaged[i] += evolution[j][i]
        time_averaged[i] /= len(evolution)
    return time_averaged

time_averaged_phi = time_averager(phi_evolution)
time_averaged_phi_stream = time_averager(phi_evolution_stream)

# Plotting
phis = np.linspace(-np.pi, np.pi, num = phi_num)
phis = np.delete(phis, 0)
phis = np.delete(phis, -1)
fig, ax = plt.subplots()
ax.plot(phis, time_averaged_phi, 
         color = 'maroon', label='Yes Stream')
ax.plot(phis, time_averaged_phi_stream,
        color = 'forestgreen', linestyle='dashed', label='No Stream')
ax.set_title(r'Time averaged $\textbf{outflow}$ density', fontsize = 18)
ax.set_xlabel(r'$\phi$ [rad] (Polar Angle)', fontsize = 14)
ax.set_ylabel(r'$\rho \left[ \frac{\mathrm{g}}{\mathrm{cm}^3} \right] $',fontsize = 14)
ax.set_yscale('log')
plt.legend(loc='upper left')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.7, 0.8, 'This uses only 3 snapshots!, \n It will be deployed at scale.',
         transform=ax.transAxes, bbox = props)
plt.grid()
    
    