#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:50:20 2023

@author: konstantinos
logρ - logT
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
AEK = '#F1C410' # Important color
import colorcet

alice = False

from src.Calculators.casters import THE_SMALL_CASTER
from src.Extractors.time_extractor import days_since_distruption
    
def pericenter_finder(radii):
    peri = 0.5
    best = 8415 # arb
    for i in range(len(radii)):
        if np.abs(radii[i] - peri) < np.abs(best - peri):
            best = radii[i]
            peri_idx = i
    return peri_idx

#%%
def maker(m, fix, pixel_num, alice):
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
    # Constants
    G = 6.6743e-11 # SI
    Msol = 1.98847e30 # kg
    Rsol = 6.957e8 # m
    t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
    # Need these for the PW potential
    c = 3e8 * t/Rsol # c in simulator units.
    rg = 2*Mbh/c**2
    # Choose snapshot
    fix = str(fix)
    
    apocenter = 2 * Rt * Mbh**(1/3)
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
    Z = np.load(folder + '/CMz_' + fix + '.npy')
    R = np.sqrt(X**2 + Y**2 + Z**2)

    # CM Position Data
    Vx = np.load(folder + '/Vx_' + fix + '.npy')
    Vy = np.load(folder + '/Vy_' + fix + '.npy')
    Vz = np.load(folder + '/Vz_' + fix + '.npy')
    V = np.sqrt(Vx**2 + Vy**2 + Vz**2)
    del X, Y, Z, Vx, Vy, Vz
    # Import Mass
    if alice:
        M = np.load(folder + '/Mass__' + fix + '.npy')
    else:
        M = np.load(folder + '/Mass_' + fix + '.npy')
    # # We only care about bound material
    Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
    bound_mask = np.where(Orbital < 0, 1, 0) # returns 1 on bound gas, 0 on unbound
    print(np.sum( np.multiply(M, bound_mask) ))
    
    # Import Density
    Den = np.load(folder + '/Den_' + fix + '.npy')
    T = np.load(folder + '/T_' + fix + '.npy')
    Msol_to_g = 1.989e33
    Rsol_to_cm = 6.957e10
    converter = Msol_to_g / Rsol_to_cm**3
    Den *=  converter
    
    Den = np.multiply(Den, bound_mask)
    T = np.multiply(T, bound_mask)
    
    days = np.round(days_since_distruption(folder + '/snap_'+fix+'.h5')/t_fall,2) 
    
    radii_start = 0.2*2*Rt
    radii_stop = apocenter
    radii = np.linspace(radii_start, radii_stop, pixel_num)
    peri_idx = pericenter_finder(radii / (2*Rt))
    
    Den_cast = THE_SMALL_CASTER(radii, R, Den, weights = M, avg = False)
    T_cast = THE_SMALL_CASTER(radii, R, T, weights = M, avg = True)
    # Orb_cast = THE_SMALL_CASTER(radii, R, Orbital, avg = True)
    
    # Fix fuckery
    Den_cast = np.log10(Den_cast) # we want a log plot
    Den_cast = np.nan_to_num(Den_cast, neginf=0)
    T_cast = np.log10(T_cast) # we want a log plot
    T_cast = np.nan_to_num(T_cast, neginf=0)
    # Orb_cast = np.nan_to_num(Orb_cast, neginf=0)

    return Den_cast, T_cast, radii/apocenter, days, peri_idx

def saver(m, fix, d, t, r, days, peri):
    path = 'products/ion/' + str(m) + '/'  + str(fix)
    np.savez(path, den=d, temp=t, radii=r, time=days, peri=peri)
#%% Doing the thing

final4 = 412
frames4 = 187
start4 = final4 - frames4

final6 = 1008
frames6 = 400
start6 = final6 - frames6

pixel_num4 = 1000
pixel_num6 = 1000

if alice:
    for i in range(frames4):
        fix4 = i + start4 + 1
        d4, t4, r4, days4, peri4 = maker(4, fix4, pixel_num4, alice)
        saver(4, fix4, d4, t4, r4, days4, peri4)
        
    for i in range(frames6):    
        fix6 = i + start6 + 1
        d6, t6, r6, days6, peri6 = maker(6, fix6, pixel_num6, alice)
        saver(6, fix6, d6, t6, r6, days6, peri6)
else:
    when = 'last' # early mid late
    if when == 'early':
        fix4 = 177
        fix6 = 683
    if when == 'mid':
        fix4 = 232 # 350 old is the good one
        fix6 = 844
    if when == 'last':
        fix4 = 263
        fix6 = 1008
    d4, t4, r4, days4, peri4 = maker(4, fix4, pixel_num4, alice)
    d6, t6, r6, days6, peri6 = maker(6, fix6, pixel_num6, alice)
    
    from src.Utilities.finished import finished
    finished()

#%% Make plots

if alice == False:
    # Fig init
    fig, ax = plt.subplots(1,2, num=1, clear=True, tight_layout = True,
                           sharex = True, sharey = True)
    
    # Image making
    img1 = ax[0].scatter(d6, t6,
                         c = r6, cmap = 'cet_bmy',
                         s=20, marker='h', zorder = 3)
    ax[0].scatter(d6[peri6], t6[peri6],
                         c = 'r', cmap = 'cet_bmy',
                         s=40, marker='x', zorder = 3)
    # plt.colorbar(img1)
    end = len(d4) # 100
    img2 = ax[1].scatter(d4[0:end], t4[0:end], 
                         c = r4[0:end], cmap = 'cet_bmy',
                         s=20, marker='h', zorder = 3)
    ax[1].scatter(d4[peri4], t4[peri4],
                         c = 'r', cmap = 'cet_bmy',
                         s=40, marker='x', zorder = 3)
    cax = fig.add_axes([0.99, 0.065, 0.02, 0.86])
    fig.colorbar(img1, cax=cax)
    
    if when == 'early':
        inset_flag = False
    else:
        inset_flag = True
    if inset_flag:
        l, b, h, w = .565, .35, .4, .145
        inset = fig.add_axes([l, b, w, h])
        inset.tick_params(axis = 'both', which = 'both', direction='in')
        instart = 300
        instop = 1_000
        instep = 1
        inset.scatter(d4[instart:instop:instep], t4[instart:instop:instep],
                    c = img2.to_rgba(r4[instart:instop:instep]),
                    s=8, marker='h', zorder = 5)
        inset.set_xlim(-6.5, -5.2)
        inset.set_ylim(4.15, 4.5)
        inset.xaxis.set_major_locator(plt.MaxNLocator(2))
        inset.yaxis.set_major_locator(plt.MaxNLocator(2))
        inset.grid()
         
    # Days text
    dayx = 0.03
    dayy = 0.03
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    txt1 = ax[0].text(dayx, dayy, 't/tfb: ' + str(days6), fontsize = 10,
    		    color='k', fontfamily = 'monospace', bbox = props,
                    transform=ax[0].transAxes)
    #txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
    txt2 = ax[1].text(dayx, dayy, 't/tfb: ' + str(days4), fontsize = 10,
    		    color='k', fontfamily = 'monospace', bbox = props,
                    transform=ax[1].transAxes)
    # # Ionized text
    # ionx = 0.69
    # iony = 0.19
    # txt1 = ax[0].text(ionx, iony, '50 \% Ionized', fontsize = 12,
    # 		    color='k', fontfamily = 'monospace', rotation = 8,
    #                 transform=ax[0].transAxes)
    # #txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
    # ionx = 0.15
    # iony = 0.15
    # txt1 = ax[1].text(ionx, iony, '50 \% Ionized', fontsize = 12,
    # 		    color='k', fontfamily = 'monospace', rotation = 8,
    #                 transform=ax[1].transAxes)
    # # Ionized text 2
    # ionx = 0.69
    # iony = 0.28
    # txt1 = ax[0].text(ionx, iony, '95 \% Ionized', fontsize = 12,
    # 		    color='k', fontfamily = 'monospace', rotation = 8,
    #                 transform=ax[0].transAxes)
    # #txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
    # ionx = 0.15
    # iony = 0.23
    # txt1 = ax[1].text(ionx, iony, '95 \% Ionized', fontsize = 12,
    # 		    color='k', fontfamily = 'monospace', rotation = 7,
    #                 transform=ax[1].transAxes)
    # #txt2.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
    
    # # Distance text 
    ionx = 1.05
    iony = 0.35
    #txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
    txt1 = fig.text(ionx, iony, 'Distance to BH [r/R$_a$]', fontsize = 12,
    		    color='k', fontfamily = 'monospace', rotation = 270)
    
    # Helium 2 99 % line
    y = 0.99
    kb = 8.617e-5 # [ev/K]
    c2 = -54.416 / kb 
    c3 = 4.01e-9  / 4 # for helium
    c4 = y**2 / (1-y)
    T_line = np.linspace(3, 12, num=200)
    
    den_line = np.log10(c3 / c4) + 1.5 * T_line + c2/10**T_line * np.log10( 2.718281) # e
    
    ax[0].plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)    
    ax[1].plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)

    ax[0].fill_between(den_line, T_line, 
                    color='lemonchiffon', zorder=1)
    ax[1].fill_between(den_line, T_line,
                    color='lemonchiffon', zorder=1)
    if inset_flag:
        inset.plot(den_line, T_line,
                 c = 'k', linestyle='dashed',
                 zorder=5)
        inset.fill_between(den_line, T_line,
                        color='lemonchiffon', zorder=1)
    
    # Helium 1 99 % line
    y = 0.99
    kb = 8.617e-5 # [ev/K]
    c2 = -24.587 / kb 
    c3 = 4.01e-9  / 4 # for helium
    c4 = y**2 / (1-y)
    T_line = np.linspace(3, 12, num=200)
    
    den_line = np.log10(c3 / c4) + 1.5 * T_line + c2/10**T_line * np.log10( 2.718281) # e
    
    ax[0].plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)    
    ax[1].plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)

    ax[0].fill_between(den_line, T_line, 
                    color='powderblue', zorder=1)
    ax[1].fill_between(den_line, T_line,
                    color='powderblue', zorder=1)
    if inset_flag:
        inset.plot(den_line, T_line,
                 c = 'k', linestyle='dashed',
                 zorder=5)
        inset.fill_between(den_line, T_line,
                        color='powderblue', zorder=1)
    
    
    # 90 % line
    c2 = -1.578e5
    c3 = 4.01e-9 # where we go from NR to ER
    y = 0.95
    c4 = y**2 / (1-y)
    T_line = np.linspace(3, 12, num=200)
    den_line = np.log10(c3 / c4) + 1.5 * T_line + c2/10**T_line * np.log10( 2.718281) # e
    
    ax[0].plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)    
    ax[1].plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)

    ax[0].fill_between(den_line, T_line, alpha = 0.2,
                    color='r', zorder=1)
    ax[1].fill_between(den_line, T_line, alpha = 0.2,
                    color='r', zorder=1)
    
    if inset_flag:
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
    
    ax[0].plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)    
    ax[1].plot(den_line, T_line,
             c = 'k', linestyle='dashed',
             zorder=5)

    ax[0].fill_between(den_line, T_line, 
                    color='#d2e7d6', zorder=1)
    ax[1].fill_between(den_line, T_line,
                    color='#d2e7d6', zorder=1)
    
    if inset_flag:
        inset.plot(den_line, T_line,
                 c = 'k', linestyle='dashed',
                 zorder=5)
        inset.fill_between(den_line, T_line,
                        color='#d2e7d6', zorder=1)
        
    # Axis labels
    fig.text(0.5, -0.01, r'$log_{10}(\rho)$  [g/cm$^3$]', ha='center', fontsize = 14)
    fig.text(-0.02, 0.5, r'$log_{10}(T)$  [K]', va='center', rotation='vertical', fontsize = 14)
    # ax[0].set_xlabel(r'$log_{10}(\rho)$  [g/cm$^3$]', fontsize = 12)
    # ax[0].set_ylabel(r'$log_{10}(T)$  [K]', fontsize = 12)
    #ax[1].set_xlabel(r'$log_{10}(\rho)$  [g/cm$^3$]', fontsize = 12)
    #ax[1].set_ylabel(r'$log_{10}(T)$  [K]', fontsize = 12)
    ax[0].tick_params(axis = 'both', which = 'both', direction='in')
    ax[1].tick_params(axis = 'both', which = 'both', direction='in')
    
    # Axis lims
    # if when == 'early':
    #     xlow = -26
    #     xhigh = -3.5
    #     ylow = -3
    #     yhigh = 9.2
    #     ax[0].set_xlim(xlow, xhigh)
    #     ax[0].set_ylim(0, yhigh)
    #     ax[1].set_xlim(xlow, xhigh)
    #     ax[1].set_ylim(ylow, yhigh)
    # else:
    xlow = -9
    xhigh = -4.5
    ylow = 3
    yhigh = 8.2
    ax[0].set_xlim(xlow, xhigh)
    ax[0].set_ylim(ylow, yhigh)
    ax[1].set_xlim(xlow, xhigh)
    ax[1].set_ylim(ylow, yhigh)
    
    # Titles
    # fig.suptitle('Ionization of debris', fontsize = 17)
    ax[0].set_title('$10^6 M_\odot$', fontsize = 15)
    ax[1].set_title('$10^4 M_\odot$', fontsize = 15)
    ax[0].grid(zorder = 3)
    ax[1].grid(zorder = 3)
        
    # Stupid lines to prove I cant use OPAL opacity tables
    OPAL = False
    if OPAL:
         logT = np.linspace(3.75, 8)
         logrho0 = -8 + 3*logT - 18 # R = -5
         logrho1 = -5 + 3*logT - 18 # R = -5
         logrho2 = 1 + 3*logT - 18 # R = -5
         ax[0].plot(logrho0, logT, c='r')
        # ax[0].plot(logrho1, logT, c='r')
         ax[0].plot(logrho2, logT, c='r')
         ax[1].plot(logrho0, logT, c='r')
         ax[1].plot(logrho2, logT, c='r')
    