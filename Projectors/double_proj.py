#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 12:28:59 2023

@author: konstantinos

side by side movies
"""

import numpy as np
import os
alice = False

if alice:
    from casters import THE_CASTER
    from time_extractor import days_since_distruption
else:
    from calculators.casters import THE_CASTER
    from extractors.time_extractor import days_since_distruption
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as PathEffects
    plt.rcParams['text.usetex'] = True
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [9 , 8]
    import colorcet

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
def maker(m, fix, pixel_num, alice):
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
    # Choose snapshot
    fix = str(fix)

    if alice:
        if m==4:
            folder = 'tde_data2/snap_' + fix
        if m==6:
            folder = 'tde_data/snap_' + fix
    else:
        folder = fix
        
    # CM Position Data
    X = np.load(folder + '/CMx_' + fix + '.npy')
    Y = np.load(folder + '/CMy_' + fix + '.npy')

    # Import Density
    Den = np.load(folder + '/Den_' + fix + '.npy')
    
    days = np.round(days_since_distruption(folder + '/snap_'+fix+'.h5')/t_fall,2) 

    # Need to convert Msol/Rsol^2 to g/cm
    Msol_to_g = 1.989e33
    Rsol_to_cm = 6.957e10
    converter = Msol_to_g / Rsol_to_cm**2
    Den *=  converter
    
    # Specify new grid:
    if m == 6:
        x_start = -100*2*Rt 
        x_stop = 4*2*Rt
        y_start = -20*2*Rt
        y_stop = 20*2*Rt
    elif m == 4:
        x_start = -21.5*2*Rt
        x_stop = 2*2*Rt
        y_start = -15*2*Rt
        y_stop = 15*2*Rt

    x_num = pixel_num # np.abs(x_start - x_stop)
    xs = np.linspace(x_start, x_stop, num = x_num )

    y_num = pixel_num # np.abs(y_start - y_stop)
    ys = np.linspace(y_start, y_stop, num = y_num) 
    
    # EVOKE
    den_cast = THE_CASTER(xs, X, ys, Y, Den)
    
    # Remove bullshit and fix things
    den_cast = np.nan_to_num(den_cast.T)
    den_cast = np.log10(den_cast) # we want a log plot
    den_cast = np.nan_to_num(den_cast, neginf=0) # fix the fuckery
    
    if m == 4:
        den_cut = 9
    elif m == 6:
        den_cut = 9
        
    # Color re-normalization
    den_cast[den_cast<0.1] = 0
    den_cast[den_cast>den_cut] = den_cut
    # den_cast = ndimage.rotate(den_cast, 0)
    return xs/(2*Rt), ys/(2*Rt), den_cast, days

def saver(m, fix, x, y, d, t):
    path = 'products/XY-Proj2/' + str(m) + '/'  + str(fix)
    np.savez(path, xs=x, ys=y, den_cast=d, time=t)

#%%
final4 = 412
frames4 = 187
start4 = final4 - frames4

final6 = 1008
frames6 = 400
start6 = final6 - frames6

pixel_num4 = 500
pixel_num6 = 1000

if alice:
    for i in range(frames4):
        fix4 = i + start4 + 1
        x4, y4, d4, t4 = maker(4, fix4, pixel_num4, alice)
        saver(4, fix4, x4, y4, d4, t4)
        
    for i in range(frames6):    
        fix6 = i + start6 + 1
        x6, y6, d6, t6 = maker(6, fix6, pixel_num6, alice)
        saver(6, fix6, x4, y4, d6, t6)
else:
    fix4 = 350
    fix6 = 820
    x4, y4, d4, t4 = maker(4, fix4, pixel_num4, alice)
    x6, y6, d6, t6 = maker(6, fix6, pixel_num6, alice)
#%%
if alice == False: 
    # Fig init
    fig, ax = plt.subplots(2,1, num=1, clear=True, tight_layout = True)
    
    # Image making
    img1 = ax[0].pcolormesh(x6, y6, d6, cmap='cet_fire')
    # plt.colorbar(img1)
    img2 = ax[1].pcolormesh(x4, y4, d4, cmap='cet_fire')
    
    cax = fig.add_axes([1, 0.045, 0.02, 0.905])
    fig.colorbar(img1, cax=cax)
    # Days text
    txt2 = ax[0].text(0.9, 0.93, 't/tfb: ' + str(t6),
    		    color='white', fontweight = 'bold', fontfamily = 'monospace',
                    transform=ax[0].transAxes)
    txt2.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='k')])
    
    txt1 = ax[1].text(0.9, 0.93, 't/tfb: ' + str(t4),
    		    color='white', fontweight = 'bold', fontfamily = 'monospace',
                    transform=ax[1].transAxes)
    txt1.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='k')])    
    
    # Axis labels
    fig.text(0.5, -0.01, r'X - Coordinate [x/2R$_T$]', ha='center', fontsize = 14)
    fig.text(-0.02, 0.5, r'Y - Coordinate [y/2R$_T$]', va='center', rotation='vertical', fontsize = 14)
    # Titles
    #fig.suptitle(r'XY $log_{10}(\rho)$ Projection $g/cm^3$', fontsize = 17)
    ax[0].set_title('$10^6 M_\odot$', fontsize = 15)
    ax[1].set_title('$10^4 M_\odot$', fontsize = 15)
    
    cbx = 1.05
    cby = 0.35
    txt1 = fig.text(cbx, cby, r'Density $\log_{10}(\rho)$ [g/cm$^3$]', fontsize = 15,
    		    color='k', fontfamily = 'monospace', rotation = 270)
    
    # # Axis lims
    # xlow = -1
    # xhigh = 1
    # ylow = -1
    # yhigh = 1
    # ax[0].set_xlim(xlow, xhigh)
    # ax[0].set_ylim(ylow, yhigh)
    # ax[1].set_xlim(xlow, xhigh)
    # ax[1].set_ylim(ylow, yhigh)

