#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 23:06:48 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import colorcet
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8 , 4]
plt.rcParams['axes.facecolor']='whitesmoke'
#%% Load

# Load ALICE data
name = 'products/'
tc4 = np.load(name + 'tc4.npy')
tc6 = np.load(name + 'tc6.npy')

# Log
tc4 = np.log10(tc4)
tc4 = np.nan_to_num(tc4, nan = 10)
tc6 = np.log10(tc6)
tc6 = np.nan_to_num(tc6)

# Load times
tfb4 = np.load('tfb4.npy')
tfb6 = np.load('tfb6.npy')
start4 = np.where(tfb4==0.5)
start4 = start4[0][0]
start6 = np.where(tfb6==0.5)
start6 = start6[0][0]
days4 = tfb4[start4:]
days6 = tfb6[start6:]

# Constants
Mbh6 = 10**6
Rt6 =  Mbh6**(1/3) # Msol = 1, Rsol = 1
apocenter6 = 2 * Rt6 * Mbh6**(1/3)
Mbh4 = 10**4
Rt4 =  Mbh4**(1/3) # Msol = 1, Rsol = 1
apocenter4 = 2 * Rt4 * Mbh4**(1/3)

# Radii
r4 = np.logspace(np.log10(0.1*2*Rt4), np.log10(apocenter4), num = 100)
r6 = np.logspace(np.log10(0.1*2*Rt6), np.log10(apocenter6), num = 100)

#%% Plotting
fig, ax = plt.subplots(1,2, tight_layout = True, sharex = True)

img = ax[0].pcolormesh(r4/(100*2*Rt4), days4, tc4,
                    cmap = 'cet_rainbow4', vmin = 4, vmax = 10)
ax[0].set_xlabel('r [r/R$_a$]')
ax[0].set_ylabel(r' Time / Fallback time $\left[ t/t_{FB} \right]$')

ax[0].set_title(r'$10^4$ $M_\odot $', fontsize = 20)

ax[1].pcolormesh(r6/(100*2*Rt6), days6, tc6, 
                    cmap = 'cet_rainbow4', vmin = 4, vmax = 10 )
ax[1].set_xlabel('r [r/R$_a$]')
ax[1].set_ylabel(r' Time / Fallback time $\left[ t/t_{FB} \right]$')
ax[1].set_title(r'$ 10^6$ $M_\odot $', fontsize = 20)

# Lims
ax[0].set_xlim(0.01, 0.2)
ax[1].set_xlim(0.01, 0.2)
ax[0].set_ylim(0.8, tfb4[-1])
ax[1].set_ylim(0.8, tfb6[-1])

# Cb txt
cbx = 1.07
cby = 0.35
txt1 = fig.text(cbx, cby, r'Cooling Time [$t_c/t_{FB}$]', fontsize = 11,
		    color='k', rotation = 270)

cax = fig.add_axes([1, 0.14, 0.035, 0.755])
fig.colorbar(img, cax=cax)