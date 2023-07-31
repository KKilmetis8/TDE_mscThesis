#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 00:37:28 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8, 4]
plt.rcParams['axes.facecolor']= 	'whitesmoke'
import colorcet

pre = 'products/cooling-time'
tc4 = np.load(pre + '/tc4.npy', allow_pickle=True)
tc6 = np.load(pre + '/tc6.npy', allow_pickle=True)
days4 = np.load(pre + '/days4.npy')
days6 = np.load(pre + '/days6.npy')
    
# Load 4
Mbh = 1e4
t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
Rt4 =  Mbh**(1/3) # Msol = 1, Rsol = 1
apocenter = 2 * Rt4 * Mbh**(1/3)

radii4_start = 0.2*2*Rt4
radii4_stop = apocenter # apocenter
radii4 = np.linspace(radii4_start, radii4_stop, 100) / apocenter 
    
# Load 6
Mbh = 1e6
t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
Rt6 =  Mbh**(1/3) # Msol = 1, Rsol = 1
apocenter = 2 * Rt6 * Mbh**(1/3)

radii6_start = 0.2*2*Rt6
radii6_stop = apocenter
radii6 = np.linspace(radii6_start, radii6_stop, 100) / apocenter
    
#%% Plotting
fig, ax = plt.subplots(1,2, tight_layout=True, sharey=True)

img1 = ax[0].pcolormesh(radii4,days4,tc4,
                cmap = 'cet_bmy', 
                norm = colors.LogNorm(vmin = 0.1, vmax = 1_000_000))
#ax[0].set_xlim(radii4[0], radii4[-5])


img2 = ax[1].pcolormesh(radii6,days6,tc6,
                cmap = 'cet_bmy',
                norm = colors.LogNorm(vmin = 0.1, vmax = 1_000_000))

cax = fig.add_axes([0.99, 0.074, 0.02, 0.835])
fig.colorbar(img1, cax=cax)

ax[0].set_ylim(0.33, 1.6)
ax[0].set_xlim(0.02, 0.2)
ax[1].set_xlim(0.02, 0.2)

#ax[0].set_xscale('log')
#ax[1].set_xscale('log')

# Distance text 
ionx = 1.06
iony = 0.3
#txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
txt1 = fig.text(ionx, iony, 'Cooling Time [$t_c/t_{FB}$]', fontsize = 14,
		    color='k', fontfamily = 'monospace', rotation = 270)

# Axis labels
fig.text(0.5, -0.01, r'r/R$_a$', ha='center', fontsize = 14)
fig.text(-0.02, 0.5, r' Time / Fallback time $\left[ t/t_{FB} \right]$', va='center', rotation='vertical', fontsize = 14)
ax[0].tick_params(axis = 'both', which = 'both', direction='in')
ax[1].tick_params(axis = 'both', which = 'both', direction='in')
ax[0].set_title(r'$10^4$ M$_\odot$')
ax[1].set_title(r'$10^6$ M$_\odot$')
# fig.suptitle('Mass Weigh Eccentricity', fontsize = 17)
#%% Check tc/t < 2
fig, ax = plt.subplots(1,2, tight_layout=True, sharey=True)

img1 = ax[0].pcolormesh(radii4,days4,tc4,
                cmap = 'cet_bmy', vmin = 0, vmax =2)
#ax[0].set_xlim(radii4[0], radii4[-5])

img2 = ax[1].pcolormesh(radii6,days6,tc6,
                cmap = 'cet_bmy',
                vmin = 0, vmax = 2)

cax = fig.add_axes([0.99, 0.074, 0.02, 0.835])
fig.colorbar(img1, cax=cax)

ax[0].set_ylim(0.33, 1.6)
ax[0].set_xlim(0.02, 0.2)
ax[1].set_xlim(0.02, 0.2)

#ax[0].set_xscale('log')
#ax[1].set_xscale('log')

# Distance text 
ionx = 1.06
iony = 0.3
#txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
txt1 = fig.text(ionx, iony, 'Cooling Time [$t_c/t_{FB}$]', fontsize = 14,
		    color='k', fontfamily = 'monospace', rotation = 270)

# Axis labels
fig.text(0.5, -0.01, r'r/R$_a$', ha='center', fontsize = 14)
fig.text(-0.02, 0.5, r' Time / Fallback time $\left[ t/t_{FB} \right]$', va='center', rotation='vertical', fontsize = 14)
ax[0].tick_params(axis = 'both', which = 'both', direction='in')
ax[1].tick_params(axis = 'both', which = 'both', direction='in')
ax[0].set_title(r'$10^4$ M$_\odot$')
ax[1].set_title(r'$10^6$ M$_\odot$')
# fig.suptitle('Mass Weigh Eccentricity', fontsize = 17)