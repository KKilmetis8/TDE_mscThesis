#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 18:19:32 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [8.0, 4.0]
import numba
from src.Eccentricity.eccentricity import  e_calc
from src.Extractors.time_extractor import days_since_distruption
from src.Calculators.casters import THE_SMALL_CASTER

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e4 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2
t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
@numba.njit
def masker(arr, mask):
    len_bound = np.sum(mask)
    new_arr = np.zeros(len_bound)
    k = 0
    for i in range(len(arr)):
        if mask[i]:
            new_arr[k] = arr[i]
            k += 1
    return new_arr
#%%
    
if __name__ == '__main__':
    fixes = ['177','232', '263'] # '820'] # '925']
    colarr = []
    fixdays = []
        
    for fix in fixes:
         X = np.load(fix + '/CMx_' + fix + '.npy')
         Y = np.load(fix + '/CMy_' + fix + '.npy')
         Z = np.load(fix + '/CMz_' + fix + '.npy')
         Vx = np.load(fix + '/Vx_' + fix + '.npy')
         Vy = np.load(fix + '/Vy_' + fix + '.npy')
         Vz = np.load(fix + '/Vz_' + fix + '.npy')
         Mass = np.load(fix + '/Mass_' + fix + '.npy')
         
         # Make Bound Mask
         R = np.sqrt( np.power(X,2) + np.power(Y,2)+ np.power(Z,2))
         V = np.sqrt( np.power(Vx,2) + np.power(Vy,2)+ np.power(Vz,2))
         Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
         bound_mask = np.where(Orbital < 0, 1, 0)

         # Apply Mask
         X = masker(X, bound_mask)
         Y = masker(Y, bound_mask)
         Z = masker(Z, bound_mask)
         # Redefine to take account only the bound
         R = np.sqrt( np.power(X,2) + np.power(Y,2)+ np.power(Z,2))
         Vx = masker(Vx, bound_mask)
         Vy = masker(Vy, bound_mask)
         Vz = masker(Vz, bound_mask)
         Bound_Mass = masker(Mass, bound_mask)
         
         position = np.array((X,Y,Z)).T # Transpose for col. vectors
         velocity = np.array((Vx,Vy,Vz)).T 
         del X, Y, Z, Vx, Vy, Vz
         
         # EVOKE eccentricity
         _ , ecc = e_calc(position, velocity, Mbh)
         
         # Cast down to 100 values
         radii = np.logspace(np.log10(17),  np.log10(2500), num = 100)
         mw_ecc_casted = THE_SMALL_CASTER(radii, R, ecc, weights = Mass,
                                        avg = True)
         # mw_ecc_casted = np.nan_to_num(mw_ecc_casted)
         colarr.append(mw_ecc_casted)
         
         day = np.round(days_since_distruption(fix +'/snap_' + fix + '.h5'),1)
         t_by_tfb = day/t_fall
         fixdays.append(t_by_tfb)
#%% Plotting
    img = plt.pcolormesh(radii,fixdays,colarr,
                   cmap = 'jet', vmin = 0, vmax = 1)
    plt.xlim(radii[0]-0.8 , radii[-12])
    plt.colorbar(img)
    plt.xscale('log')
    plt.ylabel('Days since Distruption', fontsize = 14)
    plt.xlabel(r'Radius [$R_{\odot}$]', fontsize = 14)
    plt.title('Eccentricity as a function of time and radius', fontsize = 17)
