#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 20:16:22 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
from matplotlib import ticker
import matplotlib.patheffects as PathEffects
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8, 5]
plt.rcParams['axes.facecolor']= 	'whitesmoke'
import colorcet
import numba

#%% Plotter Function
def plotter(i, m, d, t, days, peri, inset):
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    apocenter = 2 * Rt * Mbh**(1/3)
    radii_start = 0.2*2*Rt
    radii_stop = apocenter
    r = np.linspace(radii_start, radii_stop, 1000) / apocenter
    
    # Fig init
    fig, ax = plt.subplots(clear=True)

    # Image making
    img1 = ax.scatter(d, t,
                         c = r, cmap = 'cet_bmy',
                         s=20, marker='h', zorder = 3)
    # Pericentre
    ax.scatter(d[peri], t[peri],
                         c = 'r',
                         s=40, marker='x', zorder = 3)

    # cax = fig.add_axes([1, 0.065, 0.02, 0.86])
    fig.colorbar(img1)
    
    if inset:
        l, b, h, w = .565, .35, .4, .145
        inset = fig.add_axes([l, b, w, h])
        inset.tick_params(axis = 'both', which = 'both', direction='in')
        instart = 300
        instop = 1_000
        instep = 1
        inset.scatter(d[instart:instop:instep], t[instart:instop:instep],
                    c = img1.to_rgba(r[instart:instop:instep]),
                    s=8, marker='h', zorder = 5)
        inset.set_xlim(-6.5, -5.2)
        inset.set_ylim(4.2, 4.5)
        inset.xaxis.set_major_locator(plt.MaxNLocator(2))
        inset.yaxis.set_major_locator(plt.MaxNLocator(2))
        inset.grid()
    # Days text
    dayx = 0.38 
    dayy = 0.8
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    fig.text(dayx, dayy, 't/tfb: ' + str(days),
                fontsize = 10,
    		    color='k', fontfamily = 'monospace', bbox = props)
    # Ionized text 1
    # ionx = 0.69
    # iony = 0.19
    # ax.text(ionx, iony, '50 \% Ionized', fontsize = 12,
    # 		    color='k', rotation = 8,
    #                 transform=ax.transAxes)
    # # Ionized text 2
    # ionx = 0.69
    # iony = 0.28
    # ax.text(ionx, iony, '95 \% Ionized', fontsize = 12,
    # 		    color='k', rotation = 8,
    #                 transform=ax.transAxes)
    # Distance text 
    ionx = 0.85
    iony = 0.35
    fig.text(ionx, iony, 'Distance to BH [r/R$_a$]', fontsize = 12,
    		    color='k', rotation = 270)
    # 90% line
    c2 = -1.578e5
    c3 = 4.01e-9 # where we go from NR to ER
    y = 0.95
    c4 = y**2 / (1-y)
    T_line = np.linspace(3, 12, num=200)
    den_line = np.log10(c3 / c4) + 1.5 * T_line + c2/10**T_line * np.log10( 2.718281) # e
    
    ax.plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)    
    ax.fill_between(den_line, T_line, alpha = 0.2,
                    color='r', zorder=1)
    
    if inset:
        inset.plot(den_line, T_line,
                 c = 'k', linestyle='dashed',
                 zorder=5)
        inset.fill_between(den_line, T_line, alpha = 0.2,
                        color='r', zorder=1)
    # 50 % line
    y = 0.5
    c4 = y**2 / (1-y)
    T_line = np.linspace(3, 12, num=200)
    
    den_line = np.log10(c3 / c4) + 1.5 * T_line + c2/10**T_line * np.log10( 2.718281) # e
    
    ax.plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)    

    ax.fill_between(den_line, T_line, 
                    color='#d2e7d6', zorder=1)
    if inset:
        inset.plot(den_line, T_line,
                 c = 'k', linestyle='dashed',
                 zorder=5)
        inset.fill_between(den_line, T_line,
                        color='#d2e7d6', zorder=1)
        
    # Axis labels
    ax.set_xlabel(r'$log_{10}(\rho)$  [g/cm$^3$]', fontsize = 14)
    ax.set_ylabel(r'$log_{10}(T)$  [K]', fontsize = 14)
    ax.tick_params(axis = 'both', which = 'both', direction='in')

    # Axis lims
    xlow = -17
    xhigh = -1.8
    ylow = 1
    yhigh = 8.2
    ax.set_xlim(xlow, xhigh)
    ax.set_ylim(ylow, yhigh)
    
    # Titles
    if m == 4:
        ax.set_title('Ionization of debris | $10^4 M_\odot$', fontsize = 17)
    elif m==6:
        ax.set_title('Ionization of debris | $10^6 M_\odot$', fontsize = 17)
    ax.grid(zorder = 3)
    
    savepath = 'savedfigs/ion4/' + str(i)
    plt.savefig(savepath + '.png')
#%% 4
loadpath = 'products/ion'

# Read data
final4 = 412
frames4 = 187
fudge = 40
start4 = final4 - frames4
plt.ioff()
for i in range(1, frames4):
    print(np.round(i/frames4,2))
    npz4 = np.load(loadpath + '/4/' + str(start4 + i + 1) + '.npz')
    d4 = npz4['den']
    t4 = npz4['temp']
    days = npz4['time']
    peri = npz4['peri']
    plotter(i, 4, d4, t4, days, peri, inset = False)
    
#%% 6

loadpath = 'products/ion'

# Read data
final6 = 1008
frames6 = 400
start6 = final6 - frames6
plt.ioff()
for i in range(1, frames6):
    print(np.round(i/frames6,2))
    npz6 = np.load(loadpath + '/6/' + str(start6 + i + 1) + '.npz')
    d6 = npz6['den']
    t6 = npz6['temp']
    days = npz6['time']
    peri = npz6['peri']
    plotter(i, 6, d6, t6, days, peri, inset = False)