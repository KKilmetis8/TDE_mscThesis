#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 02:37:18 2023

@author: konstantinos
"""

import numpy as np
import gc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import colorcet
# Pretty plots
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [8 , 4]

m = 4 # 6 4
choice = 'YZ' # XY XZ YZ

# Constants
Mbh = 10**m
Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
apocenter = 2 * Rt * Mbh**(1/3)
if m == 6:
    fixes = ['683' , '844', '1008']
    days =  ['0.5', '1' , '1.6']  
    name = 'products/final-proj/' + choice + '-' + str(m) + '-'

    if choice == 'XY':
        x_start = -15_000
        x_stop = 1_500
        x_num = 2000 # np.abs(x_start - x_stop)
        xs = np.linspace(x_start, x_stop, num = x_num )
        y_start = -4_000
        y_stop = 4_000
        y_num = 2000 # np.abs(y_start - y_stop)
        ys = np.linspace(y_start, y_stop, num = y_num)
    if choice == 'XZ':
        x_start = -1_500
        x_stop = 1_000
        x_num = 1000 # np.abs(x_start - x_stop)
        xs = np.linspace(x_start, x_stop, num = x_num )
        y_start = -1_000
        y_stop = 1_000
        y_num = 1000 # np.abs(y_start - y_stop)
        ys = np.linspace(y_start, y_stop, num = y_num)
    if choice == 'YZ':
        x_start = -2_500
        x_stop = 2_000
        x_num = 2000 # np.abs(x_start - x_stop)
        xs = np.linspace(x_start, x_stop, num = x_num )
        y_start = -1_000
        y_stop = 1_000
        y_num = 2000 # np.abs(y_start - y_stop)
        ys = np.linspace(y_start, y_stop, num = y_num)
if m==4:
    fixes = ['177', '232', '263']
    days =  ['0.5', '1', '1.3']  
    name = 'products/final-proj/' + choice + '-' + str(m) + '-'
    if choice == 'XY':
        # Specify new grid:
        x_start =  -apocenter - 4 *2*Rt
        x_stop =  12 * 2*Rt
        x_num = 500 # np.abs(x_start - x_stop)
        xs = np.linspace(x_start, x_stop, num = x_num )
        # y +- 150, z +- 50
        y_start = -apocenter
        y_stop =  apocenter
        y_num = 500 # np.abs(y_start - y_stop)
        ys = np.linspace(y_start, y_stop, num = y_num )
    if choice == 'YZ':
        # y +- 150, z +- 50
        y_start = -apocenter
        y_stop = apocenter
        y_num = 600 # np.abs(y_start - y_stop)
        xs = np.linspace(y_start, y_stop, num = y_num)    
        z_start = -apocenter
        z_stop = apocenter
        z_num = 600 # np.abs(y_start - y_stop)
        ys = np.linspace(z_start, z_stop, num = z_num)    

for i in range(len(fixes)):
    
    den_cast = np.load(name + fixes[i] + '.npy')
    # Plotting
    fig, ax = plt.subplots()
    img = ax.pcolormesh(xs, ys, den_cast, cmap='cet_fire',
                        vmin = 0, vmax = 8)
    fig.colorbar(img)
    
    # ax.set_ylim(-50, 50)
    # ax.set_xlim(-1500, 1000)
    ax.set_title( choice + r' Density $\log_{10}( \rho )$ [g/cm$^2$] projection',
                 fontsize = 20)
    
    ax.set_xlabel(choice[0]+'-coordinate $[R_{\odot}$]')
    ax.set_ylabel(choice[1]+'-coordinate $[R_{\odot}$]')

    txt_x = xs[0] + 50
    txt_y = ys[0] + 50
    
    ax.text(txt_x, txt_y, 'Days after disruption: ' + days[i],
            color='white', 
            fontweight = 'bold', 
            fontname = 'Consolas',
            fontsize = 18)
