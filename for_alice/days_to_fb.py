#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 01:07:49 2023

@author: konstantinos
"""

import numpy as np

from time_extractor import days_since_distruption


def maker(m, fix):
    # BH specific constants
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
    # rg = 2*Mbh/c**2
    
    # Choose snapshot
    fix = str(fix)

    if m==4:
        folder = 'tde_data2/snap_' + fix
    if m==6:
        folder = 'tde_data/snap_' + fix
        
    days = np.round(days_since_distruption(folder + '/snap_'+fix+'.h5')/t_fall,2) 
    
    return days
#%% Get Z, rho, T
final4 = 412
frames4 = 187
start4 = final4 - frames4

final6 = 1008
frames6 = 400
start6 = final6 - frames6

pixel_num = 100
max_z6 = 200
max_z4 = 50

days4_store = []
days6_store = []
for i in range(frames4):
    fix4 = i + start4 + 1
    # Cast Down
    days4 = maker(4, fix4)
    days4_store.append(days4)
    
for i in range(frames6):    
    fix6 = i + start6 + 1
    # Cast Down
    days6 = maker(6, fix6)

    days6_store.append(days6)
    
# Save
savepath = 'products/cooling-time/'

np.save(savepath + 'days4', days4_store)
np.save(savepath + 'days6', days6_store)