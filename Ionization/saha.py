#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 16:53:53 2023

@author: konstantinos
"""

import numpy as np
import os
alice = False

if alice:
    from src.Calculators.casters import THE_CASTER
    from src.Extractors.time_extractor import days_since_distruption
else:
    from src.Calculators.casters import THE_CASTER
    from src.Extractors.time_extractor import days_since_distruption
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

me = 9.1e-28 # cgs
hbar =  6.662e-27 # # cgs
kb = 1.38e-16 # cgs

# Ionization energies
ev_to_erg = 1.602e-12
xh = 13.598 * ev_to_erg # erg
prefactor_h = 1 # 2g1/g0

xhe1 = 24.587 * ev_to_erg # erg
prefactor_he1 = 4

xhe2 = 54.416 * ev_to_erg # erg
prefactor_he2 = 1
 
#%%
def maker(m, fix, pixel_num, alice):
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
    apocenter = 2 * Rt * Mbh**(1/3)
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

    # Import Pressure, Tempreature
    P = np.load(folder + '/P_' + fix + '.npy')
    T = np.load(folder + '/T_' + fix + '.npy')
    
    days = np.round(days_since_distruption(folder + '/snap_'+fix+'.h5')/t_fall,2) 

    # Need to convert pressure to cgs
    Msol_to_g = 1.989e33
    Rsol_to_cm = 6.957e10
    converter = Msol_to_g / (Rsol_to_cm * t**2)
    P *=  converter

    # Caclculate Ks
    # NOTE: Add degeneracy factors
    K1 = prefactor_h * (2*np.pi/me)**1.5 * hbar**3 * np.exp(xh/(kb * T)) / (kb * T**2.5)
    ion1 = np.divide(1, np.sqrt(1 + P*K1))
    
    K2 = prefactor_he1 * (2*np.pi/me)**1.5 * hbar**3 * np.exp(xhe1/(kb * T)) / (kb * T**2.5)
    ion2 = np.divide(1, np.sqrt(1 + P*K2))
    
    K3 = prefactor_he2 * (2*np.pi/me)**1.5 * hbar**3 * np.exp(xhe2/(kb * T)) / (kb * T**2.5)
    ion3 = np.divide(1, np.sqrt(1 + P*K3))
    del K1, K2, K3, P, T
    # Ion fractions
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
    ion1_cast = THE_CASTER(xs, X, ys, Y, ion1, avg = True)
    ion1_cast = np.nan_to_num(ion1_cast.T)
    
    ion2_cast = THE_CASTER(xs, X, ys, Y, ion2, avg = True)
    ion2_cast = np.nan_to_num(ion2_cast.T)
    
    ion3_cast = THE_CASTER(xs, X, ys, Y, ion3, avg = True)
    ion3_cast = np.nan_to_num(ion3_cast.T)


    return xs/(2*Rt), ys/(2*Rt), ion1_cast, ion2_cast, ion3_cast, days

def saver(m, fix, x, y, d, t):
    path = 'products/ion_fraction/' + str(m) + '/'  + str(fix)
    np.savez(path, xs=x, ys=y, den_cast=d, time=t)

#%%
final4 = 412
frames4 = 187
start4 = final4 - frames4

final6 = 1008
frames6 = 400
start6 = final6 - frames6

pixel_num4 = 500
pixel_num6 = 700

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
    when = 'last' # early mid late
    if when == 'early':
        fix4 = 177
        fix6 = 683
    if when == 'mid':
        fix4 = 232
        fix6 = 844
    if when == 'last':
        fix4 = 263
        fix6 = 1008

    x4, y4, h4, he14, he24, t4 = maker(4, fix4, pixel_num4, alice)
    x6, y6, h6, he16, he26, t6 = maker(6, fix6, pixel_num6, alice)
    import src.Utilities.finished
#%%    
if alice == False:     
    # Fig init
    fig, ax = plt.subplots(2,1, num=1, clear=True, tight_layout = True)
    
    # Image making
    img1 = ax[0].pcolormesh(x6, y6, h6, cmap='cet_fire')
    # plt.colorbar(img1)
    img2 = ax[1].pcolormesh(x4, y4, h4, cmap='cet_fire')
    
    cax = fig.add_axes([1, 0.045, 0.02, 0.86])
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
    fig.suptitle(r'Hydrogen Ionization', fontsize = 17)
    ax[0].set_title('$10^6 M_\odot$', fontsize = 15)
    ax[1].set_title('$10^4 M_\odot$', fontsize = 15)
    
    cbx = 1.07
    cby = 0.35
    txt1 = fig.text(cbx, cby, 'Ionization Fraction', fontsize = 15,
    		    color='k', fontfamily = 'monospace', rotation = 270)
#%% ################# Helium 1st Ionization plot ##############################
    plt.figure()
    fig, ax = plt.subplots(2,1, num=1, clear=True, tight_layout = True)
    
    # Image making
    img1 = ax[0].pcolormesh(x6, y6, he16, cmap='cet_fire')
    # plt.colorbar(img1)
    img2 = ax[1].pcolormesh(x4, y4, he14, cmap='cet_fire')
    
    cax = fig.add_axes([0.93, 0.123, 0.04, 0.76])
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
    fig.text(0.5, 0.06, r'X - Coordinate [x/2R$_T$]', ha='center', fontsize = 14)
    fig.text(0.05, 0.5, r'Y - Coordinate [y/2R$_T$]', va='center', rotation='vertical', fontsize = 14)
    # Titles
    fig.suptitle(r'Helium 1st Ionization', fontsize = 17)
    ax[0].set_title('$10^6 M_\odot$', fontsize = 15)
    ax[1].set_title('$10^4 M_\odot$', fontsize = 15)
    
    cbx = 1.02
    cby = 0.45
    txt1 = fig.text(cbx, cby, 'Ionization Fraction', fontsize = 15,
    		    color='k', fontfamily = 'monospace', rotation = 270)
#%% ################# Helium 2nd Ionization plot ##############################
    plt.figure()
    fig, ax = plt.subplots(2,1, num=1, clear=True, tight_layout = True)
    
    # Image making
    img1 = ax[0].pcolormesh(x6, y6, he26, cmap='cet_fire')
    # plt.colorbar(img1)
    img2 = ax[1].pcolormesh(x4, y4, he24, cmap='cet_fire')
    
    cax = fig.add_axes([0.93, 0.123, 0.04, 0.76])
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
    fig.text(0.5, 0.06, r'X - Coordinate [x/2R$_T$]', ha='center', fontsize = 14)
    fig.text(0.05, 0.5, r'Y - Coordinate [y/2R$_T$]', va='center', rotation='vertical', fontsize = 14)
    # Titles
    fig.suptitle(r'Helium 2nd Ionization', fontsize = 17)
    ax[0].set_title('$10^6 M_\odot$', fontsize = 15)
    ax[1].set_title('$10^4 M_\odot$', fontsize = 15)
    
    cbx = 1.02
    cby = 0.45
    txt1 = fig.text(cbx, cby, 'Ionization Fraction', fontsize = 15,
    		    color='k', fontfamily = 'monospace', rotation = 270)