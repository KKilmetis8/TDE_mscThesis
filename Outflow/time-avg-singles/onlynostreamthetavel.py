#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 13:22:32 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
import gc
# Pretty Plots
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8 , 4]
plt.rcParams['axes.facecolor']='whitesmoke'

# Constants
Rsol_to_cm = 6.957e10

# Coord choices r-, theta-, phi-
coord = 'theta-'
# den, vel
quantity = 'vel'
# Dont show plots
plt.ioff()

# Data Load
intervals =  ['36-45', '45-55', '55+']
stream = ['-yes', '-no']
days = np.load('tavg-data/days.npy')
if coord=='theta-':
    # Theta range
    theta_num = 100
    thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
    thetas = np.delete(thetas, 0)
    thetas = np.delete(thetas, -1)
# θ - Velocity --------------------------------------------------------------------------                     
    if quantity == 'vel':
        thing = coord+quantity
        k = 0 # Counter for plot names
        for interval in intervals:
            ta_rden = np.load('tavg-data/evo'+thing+stream[0]+interval + '.npy')
            # ta_rden_stream = np.load('tavg-data/evo'+thing+stream[1]+interval +'.npy')
            
            # Generate thetas
            theta_num = len(ta_rden[0]) + 2
            thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
            thetas = np.delete(thetas, 0)
            thetas = np.delete(thetas, -1)
            for i in range(len(ta_rden)):
                # Plot
                fig, ax = plt.subplots(1, clear=True)
                ax.plot(thetas, ta_rden[i], 
                          color = 'navy', label='All outflows')
                # ax.plot(thetas, ta_rden_stream[i], 
                #          color = 'darkorange', linestyle='--', label='Excluding the stream')
                ax.set_title(r'$\textbf{Outflow}$ velocity for days 21+', fontsize = 18)
                ax.set_xlabel(r'Latitudinal Angle $\left[ \theta \right]$ ', fontsize = 14)
                ax.set_ylabel(r'$|V| \left[ \frac{\mathrm{km}}{\mathrm{s}} \right] $',fontsize = 14)
                ax.set_yscale('log')
                plt.legend(loc='upper center')
                plt.grid()
                ax.tick_params(axis = 'both', which = 'both', direction='in')

                # Axes Limit
                ax.set_ylim(1e0,6e4)
                
                # # Light speed
                # # ax.axhline(300_000 , c='k')
                # # plt.text(-1.5, 300_000, 'V = c', # very nice box
                # #           bbox = props)
                
                # # Day counter
                text = 'Day: ' + str(days[k])
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                plt.text(0.45, 0.1, text, # very nice box
                          transform=ax.transAxes, bbox = props)
                
                # Save
                plt.savefig('savedfigs/evo/' + thing + str(k)+ '.png')
                k +=1 