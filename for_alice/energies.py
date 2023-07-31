#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 17:44:39 2023

@author: konstantinos
"""

import numpy as np
from time_extractor import days_since_distruption

# Choose sims
bh = '6' # or '4
if bh == '6':
    Mbh = 1e6 # * Msol
    pre = 'tde_data/snap_'
    fixes = np.arange(750,1008+1) # 750,1008+1 , 225,412+1
if bh == '4':
    Mbh = 1e4 # * Msol
    pre = 'tde_data2/snap_'
    fixes = np.arange(225,412+1) # 750,1008+1 , 225,412+1
    
# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1

# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2

# Fallback rate
t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13

#%% Do the thing

orbs = []
intes = []
rads = []
times = []
for fix in fixes:
     fix = str(fix)
     X = np.load(pre + fix + '/CMx_' + fix + '.npy')
     Y = np.load(pre + fix + '/CMy_' + fix + '.npy')
     Z = np.load(pre + fix + '/CMz_' + fix + '.npy')
     Vx = np.load(pre + fix + '/Vx_' + fix + '.npy')
     Vy = np.load(pre + fix + '/Vy_' + fix + '.npy')
     Vz = np.load(pre + fix + '/Vz_' + fix + '.npy')
     Internal = np.load(pre + fix + '/IE_' + fix + '.npy')
     Rad_per_m = np.load(pre + fix + '/Rad_' + fix + '.npy')
     Mass = np.load(pre + fix + '/Mass__' + fix + '.npy')
     
     # Get time
     day = np.round(days_since_distruption(pre + fix +'/snap_' + fix + '.h5'),1)
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
# Save
energies = [times, orbs, intes, rads]
np.save('Energies-' + bh, energies)
