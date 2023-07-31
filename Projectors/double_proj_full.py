#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:34:15 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patheffects as PathEffects
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [9 , 8]
import colorcet as cc # cooler colormaps

loadpath = 'products/XY-Proj2/'
# Read data
final4 = 412
frames4 = 187
fudge = 40
start4 = final4 - frames4 + fudge

final6 = 1008
frames6 = 400
start6 = final6 - frames6

# Re-define grid because ALICE was drunk last night
m = 6
Mbh = 10**m
Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
x_start = -100*2*Rt 
x_stop = 4*2*Rt
y_start = -20*2*Rt
y_stop = 20*2*Rt
x6 = np.linspace(x_start, x_stop, 1000) / (2*Rt)
y6 = np.linspace(y_start, y_stop, 1000) / (2*Rt)

m = 4
Mbh = 10**m
Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
x_start = -21.5*2*Rt
x_stop = 2*2*Rt
y_start = -15*2*Rt
y_stop = 15*2*Rt 
x4 = np.linspace(x_start, x_stop, 500) / (2*Rt)
y4 = np.linspace(y_start, y_stop, 500) / (2*Rt)
#%% Model Loops

# 4 lags
# counter = 1
# rate = (frames4 - fudge) / frames6
# j = 0
# for i in range(1, frames6 + 1):
#     print(start6 + i)
#     if counter == 1:
#         print('-----' , start4 + j)
#         j += 1
#         counter -= 1
#     elif counter>1:
#         print('-----' , start4 + j)
#         j += 1
#         counter -= 1
#     counter += rate

# 6 skips
counter = 1
rate = frames6 / (frames4 - fudge)
j = 1
for i in range(1, frames4 - fudge + 1):
    print(start4 + i)

    print('---', start6 + int(j))
    j += rate
#%% 6 skips
plt.ioff()
counter = 1
rate = frames6 / (frames4 - fudge)
j = 1
for i in range(1, frames4 - fudge + 1):
    print(start4 + i)

    print('---', start6 + int(j))
    
    fix4 = '4/' + str(start4 + i) 
    npz4 = np.load(loadpath + fix4 + '.npz')
    # x4 = npz4['xs']
    # y4 = npz4['ys']
    d4 = npz4['den_cast']
    t4 = npz4['time']
    
    fix6 = '6/' + str(start6 + int(j)) 
    npz6 = np.load(loadpath + fix6 + '.npz')
    # x6 = npz6['xs']
    # y6 = npz6['ys']
    d6 = npz6['den_cast']
    t6 = npz6['time']

    j += rate
    # Fig init
    fig, ax = plt.subplots(2,1, num=1, clear=True, tight_layout = True)
    
    # Image making
    img1 = ax[0].pcolormesh(x6, y6, d6, cmap='cet_fire', vmax = 8, vmin = 0)
    # plt.colorbar(img1)
    img2 = ax[1].pcolormesh(x4, y4, d4, cmap='cet_fire', vmax = 8, vmin = 0)
    
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
    
    savepath = 'savedfigs/double4/' + str(i)
    plt.savefig(savepath + '.png', bbox_inches = 'tight')
#%% 4 lags
plt.ioff()
counter = 1
rate = (frames4 - fudge) / frames6
j = 0
for i in range(1): #1, frames6 + 1):
    print(start6 + i)
    if counter == 1:
        print('-----' , start4 + j)
        j += 1
        counter -= 1
        fix4 = '4/' + str(start4 + j) 
        npz4 = np.load(loadpath + fix4 + '.npz')
        # x4 = npz4['xs']
        # y4 = npz4['ys']
        d4 = npz4['den_cast']
        t4 = npz4['time']
    elif counter>1:
        print('-----' , start4 + j)
        j += 1
        counter -= 1
        fix4 = '4/' + str(start4 + j) 
        npz4 = np.load(loadpath + fix4 + '.npz')
        # x4 = npz4['xs']
        # y4 = npz4['ys']
        d4 = npz4['den_cast']
        t4 = npz4['time']    
    counter += rate
    
    fix6 = '6/' + str(start6 + i + 1) 
    npz6 = np.load(loadpath + fix6 + '.npz')
    # x6 = npz6['xs']
    # y6 = npz6['ys']
    d6 = npz6['den_cast']
    t6 = npz6['time']

    # Fig init
    fig, ax = plt.subplots(2,1, num=1, clear=True, tight_layout = True)
    
    # Image making
    img1 = ax[0].pcolormesh(x6, y6, d6, cmap='cet_fire', vmax = 8, vmin = 0)
    # plt.colorbar(img1)
    img2 = ax[1].pcolormesh(x4, y4, d4, cmap='cet_fire', vmax = 8, vmin = 0)
    
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
    
    savepath = 'savedfigs/double4/' + str(i)
    plt.savefig(savepath + '.png', bbox_inches = 'tight')