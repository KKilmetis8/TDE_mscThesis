#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 13:24:54 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [8.0, 4.0]
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.facecolor']='whitesmoke'
from extractors.time_extractor import days_since_distruption

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2

t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13

#%%
    
fixes = ['820', '881']
orbs = []
intes = []
rads = []
times = []
for fix in fixes:
     X = np.load(fix + '/CMx_' + fix + '.npy')
     Y = np.load(fix + '/CMy_' + fix + '.npy')
     Z = np.load(fix + '/CMz_' + fix + '.npy')
     Vx = np.load(fix + '/Vx_' + fix + '.npy')
     Vy = np.load(fix + '/Vy_' + fix + '.npy')
     Vz = np.load(fix + '/Vz_' + fix + '.npy')
     Internal = np.load(fix + '/IE_' + fix + '.npy')
     Rad_per_m = np.load(fix + '/Rad_' + fix + '.npy')
     Mass = np.load(fix + '/Mass_' + fix + '.npy')
     
     # Get time
     day = np.round(days_since_distruption(fix+'/snap_'+fix+'.h5'),1)
     t_by_tfb = day/t_fall
     times.append(t_by_tfb)
     
     # Calc orbital energy 
     R = np.sqrt( np.power(X,2) + np.power(Y,2)+ np.power(Z,2))
     V = np.sqrt( np.power(Vx,2) + np.power(Vy,2)+ np.power(Vz,2))
     Orbital = (0.5 * V**2 ) + Mbh / (R-rg) # positive
     del X, Y, Z, Vx, Vy, Vz
     
     # Calc rad renergy
     Rad = np.multiply(Mass, Rad_per_m) 
     
     # Cast down to 100 values
     orb = Orbital.sum()
     inte = Internal.sum()
     rad = Rad.sum()
     
     orbs.append(orb)
     intes.append(inte)
     rads.append(rad)

#%% Plotting

plt.plot(times, orbs,  '-o', c='tab:purple', label='Orbital')
plt.plot(times, intes, '-s', c='tab:red', label = 'Internal')
plt.plot(times, rads, '-^', c='tab:olive', label = 'Radiation')

plt.yscale('log')
plt.tick_params(axis = 'both', which = 'both', direction='in')
plt.legend()
plt.grid()
plt.ylabel(r'Total Energy [$M_{\odot} R_{\odot}^2 \widetilde{t}^{-2}$]', fontsize = 14)
plt.xlabel(r' Time / Fallback time $\left[ t/t_{FB} \right]$', fontsize = 14)
plt.title('Energy evolution', fontsize = 17)
