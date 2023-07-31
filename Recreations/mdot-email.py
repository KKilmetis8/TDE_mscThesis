#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 08:48:42 2023

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
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
Rsol_to_cm = 6.957e10

converter = Rsol**2 * 1e6 # kg/s [SI] See: notes p.33
to_msol = converter/Msol

# Coord choices r-, theta-, phi-
coord = 'r-'
# den, vel
den = 'den'
vel = 'vel'
thing1 = coord + den
thing2 = coord + vel
# Dont show plots
plt.ioff()

# Data Load
intervals =  ['36-45', '45-55', '55+']
stream = ['-yes', '-no']
days = np.load('tavg-data/days700.npy')

# Radii range 
# NEEDS TO BE THE SAME ONE YOU USED TO MAKE THE CASTINGS
r_start = 50
r_stop = 80_000
r_num = 100
radii = np.linspace(r_start, r_stop, num = r_num)
# First and last points always get smooshed.
radii = np.delete(radii, 0) 
radii = np.delete(radii, -1)

#--- The loop begins --------------------------------------------------------------------------                             
k = 0 # Counter for plot names

for interval in intervals:
    # Density
    ta_rden = np.load('tavg-data/evo' + thing1 + stream[0] + interval + '.npy')
    ta_rden_stream = np.load('tavg-data/evo' + thing1 + stream[1] + interval +'.npy')
    # Velocity
    ta_rvel = np.load('tavg-data/evo' + thing2 + stream[0] + interval + '.npy')
    ta_rvel_stream = np.load('tavg-data/evo' + thing2 + stream[1] + interval +'.npy')
    for i in range(2): #len(ta_rden)):
        # Calc Mdot = 4πr^2 V(r) ρ(r)
        # mdot = 4 * np.pi * radii**2 * ta_rden[i] * ta_rvel[i]
        mdot_stream = 4*np.pi*radii**2 * ta_rden_stream[i] * ta_rvel_stream[i]
        
        # # Convert to Msol/s
        # mdot *= to_msol
        mdot_stream *= to_msol
        
        fig, ax = plt.subplots(1, clear=True)
        # ax.plot(radii, mdot, 
        #           color = 'purple', label='All outflows')
        ax.plot(radii, mdot_stream, 
                color = 'goldenrod', linestyle='--', label='Excluding the stream')
        ax.set_title(r'$\textbf{Outflowing}$ $\dot{M}$ for days 21.5+', fontsize = 18)
        ax.set_xlabel(r'Radial Distance $\left[ R_{\odot} \right]$ ', fontsize = 14)
        ax.set_ylabel(r'$\dot{M} \left[ \frac{\mathrm{M_{\odot}}}{\mathrm{s}} \right] $',fontsize = 14)
        ax.set_yscale('log')
        plt.legend(loc='upper right')
        plt.grid()
        ax.tick_params(axis = 'both', which = 'both', direction='in')
        
        # # Axis limit
        plt.ylim(1e-10, 1e0)

        # # Day counter
        text = 'Day: ' + str(days[k])
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        plt.text(0.0465, 0.9, text, # very nice box
                  transform=ax.transAxes, bbox = props)
        
        # Save
        plt.savefig('savedfigs/evo/' + 'mdot' + str(k)+ '.png')
        k +=1 
        gc.collect()
                
