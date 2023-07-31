#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 15:08:52 2023

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

# Coord choices r-, theta-, phi-
m = '4'
coord = 'r-'
# den, vel
quantity = 'vel'
Mbh = '-M' + m + coord
# Dont show plots
plt.ioff()
stream = ['no']

if m=='4':
    Mbh = 1e4
    t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    apocenter = 2 * Rt * Mbh**(1/3)
    tfb = np.load('tfb4.npy')
if m=='6':
    Mbh = 1e6
    t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    apocenter = 2 * Rt * Mbh**(1/3)
    tfb = np.load('tfb6.npy')

days = np.load('tavg-data/4/days4.npy')
if coord == 'r-':
    # Radii range 
    # NEEDS TO BE THE SAME ONE YOU USED TO MAKE THE CASTINGS
    r_start = 0.5*Rt
    r_stop = apocenter
    r_num = 200
    radii = np.linspace(r_start, r_stop, num = r_num)

# r - Density --------------------------------------------------------------------------                             
    if quantity == 'den':
        k = 0 # Counter for plot names
        thing = coord+quantity
        for part in range(2):
            part = '-' + str(part+1)
            ta_rden = np.load('tavg-data/evo'+Mbh+thing+stream[0]+part + '.npy')
            ta_rden_stream = np.load('tavg-data/evo'+Mbh+thing+stream[1]+part +'.npy')

            for i in range(len(ta_rden)):
                # Plot
                fig, ax = plt.subplots(1, clear=True)
                ax.plot(radii, ta_rden[i], 
                         color = 'purple', label='All outflows')
                ax.plot(radii, ta_rden_stream[i], 
                         color = 'goldenrod', linestyle='--', label='Excluding the stream')
                ax.set_title(r'Mass-Weighted $\textbf{Outflow}$ density for days 21.5+', fontsize = 18)
                ax.set_xlabel(r'Radial Distance $\left[ R_{\odot} \right]$ ', fontsize = 14)
                ax.set_ylabel(r'$\rho \left[ \frac{\mathrm{g}}{\mathrm{cm}^3} \right] $',fontsize = 14)
                ax.set_yscale('log')
                plt.legend(loc='upper right')
                plt.grid()
                ax.tick_params(axis = 'both', which = 'both', direction='in')
                
                # Axis limit
                # plt.ylim(1e-22, 1e-7)

                # Day counter
                text = 'Day: ' + str(days[k])
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                plt.text(0.0465, 0.9, text, # very nice box
                         transform=ax.transAxes, bbox = props)
                
                # Save
                plt.savefig('savedfigs/evo/' + thing + str(k)+ '.png')
                k +=1 
                gc.collect()
                
# r - Velocity --------------------------------------------------------------------------                         
    if quantity == 'vel':
        thing = Mbh+quantity
        k = 0 # Counter for plot names
        for part in range(2):
            part = '-' + str(part+1)
            ta_rden = np.load('tavg-data/evo'+thing+stream[0]+part + '.npy')
            # ta_rden_stream = np.load('tavg-data/evo'+thing+stream[1]+part +'.npy')
            for i in range(len(ta_rden)):
                # Plot
                fig, ax = plt.subplots(1, clear=True)
                ax.plot(radii, ta_rden[i], 
                         color = 'purple', label='All outflows')
                # ax.plot(radii, ta_rden_stream[i], 
                #          color = 'goldenrod', linestyle='--', label='Excluding the stream')
                ax.set_title(r'Mass-Weighted $\textbf{Outflow}$ velocity for days 21.5+', fontsize = 18)
                ax.set_xlabel(r'Radial Distance $\left[ R_{\odot} \right]$ ', fontsize = 14)
                ax.set_ylabel(r'$|V| \left[ \frac{\mathrm{km}}{\mathrm{s}} \right] $',fontsize = 14)
                ax.set_yscale('log')
                plt.legend(loc='lower right')
                ax.tick_params(axis = 'both', which = 'both', direction='in')
                plt.grid()
                
                # Axes Limit
                ax.set_ylim(1e-2,1e5)
                
                # Day counter
                text = 'Day: ' + str(days[k])
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                plt.text(0.0465, 0.9, text, # very nice box
                         transform=ax.transAxes, bbox = props)
                
                # Save
                plt.savefig('savedfigs/evo/' + thing + str(k)+ '.png')
                k +=1 
if coord=='theta-':
    # Theta range
    theta_num = 100
    thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
# θ - Density --------------------------------------------------------------------------                         
    if quantity == 'den':
        thing = coord+quantity
        k = 0 # Counter for plot names
        for part in range(4):
            part = '-' + str(part+1)
            ta_rden = np.load('tavg-data2/evo'+thing+stream[0]+part+'.npy')
            ta_rden_stream = np.load('tavg-data2/evo'+thing+stream[1]+part+'.npy')            
            for i in range(2): #len(ta_rden)):
                # Plot
                fig, ax = plt.subplots(1, clear=True)
                ax.plot(thetas, ta_rden[i], 
                         color = 'navy', label='All outflows')
                ax.plot(thetas, ta_rden_stream[i], 
                         color = 'darkorange', linestyle='--', label='Excluding the stream')
                ax.set_title(r'Mass-Weighted $\textbf{Outflow}$ density for days 21+', fontsize = 18)
                ax.set_xlabel(r'Latitudinal Angle $\left[ \theta \right]$ ', fontsize = 14)
                ax.set_ylabel(r'$\rho \left[ \frac{\mathrm{g}}{\mathrm{cm}^3} \right] $',fontsize = 14)
                ax.set_yscale('log')
                plt.legend(loc='upper right')
                plt.grid()
                ax.tick_params(axis = 'both', which = 'both', direction='in')      
                
                # Axes Limit
                ax.set_ylim(1-15,2e-8)
                
                # Day Ticker
                text = 'Day: ' + str(days[k])
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                plt.text(0.0465, 0.9, text, # very nice box
                         transform=ax.transAxes, bbox = props)
                
                # Save
                plt.savefig('savedfigs/evo/'+thing+str(k)+'.png')
                k += 1
                
                gc.collect()
# θ - Velocity --------------------------------------------------------------------------                     
    if quantity == 'vel':
        thing = coord+quantity
        k = 0 # Counter for plot names
        for part in range(4):
            part = '-' + str(part+1)
            ta_rden = np.load('tavg-data2/evo'+thing+stream[0]+part+'.npy')
            ta_rden_stream = np.load('tavg-data2/evo'+thing+stream[1]+part+'.npy')
            
            # Generate thetas
            theta_num = len(ta_rden[0])
            thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
            for i in range(2): #len(ta_rden)):
                # Plot
                fig, ax = plt.subplots(1, clear=True)
                ax.plot(thetas, ta_rden[i], 
                          color = 'navy', label='All outflows')
                ax.plot(thetas, ta_rden_stream[i], 
                          color = 'darkorange', linestyle='--', label='Excluding the stream')
                # Surroundings
                ax.set_title(r'Mass-Weighted $\textbf{Outflow}$ velocity for days 21+', fontsize = 18)
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
                
                # Day counter
                text = 'Day: ' + str(days[k])
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                plt.text(0.45, 0.1, text, # very nice box
                          transform=ax.transAxes, bbox = props)
                
                # Save
                plt.savefig('savedfigs/evo/' + thing + str(k)+ '.png')
                k +=1 
if coord=='phi-':
    # Phi range
    phi_num = 100
    phis = np.linspace(-np.pi, np.pi, num = phi_num)
# φ - Density --------------------------------------------------------------------------     
    if quantity == 'den':
        thing = coord+quantity
        k = 0
        for part in range(4):
            part = '-' + str(part+1)
            ta_rden = np.load('tavg-data2/evo'+thing+stream[0]+part+ '.npy')
            ta_rden_stream = np.load('tavg-data2/evo'+thing+stream[1]+part+'.npy')
            for i in range(len(ta_rden)):
                # Plot
                fig, ax = plt.subplots(1, clear=True)
                ax.plot(phis, ta_rden[i], 
                         color = 'maroon', label='All outflows')
                ax.plot(phis, ta_rden_stream[i], 
                         color = 'forestgreen', linestyle='--', label='Excluding the stream')
                ax.set_title(r'Mass-Weighted $\textbf{Outflow}$ density for days 21+', fontsize = 18)
                ax.set_xlabel(r'Longitudinal Angle $\left[ \phi \right]$ ', fontsize = 14)
                ax.set_ylabel(r'$\rho \left[ \frac{\mathrm{g}}{\mathrm{cm}^3} \right] $',fontsize = 14)
                ax.set_yscale('log')
                plt.legend(loc='upper right')
                ax.tick_params(axis = 'both', which = 'both', direction='in')   
                plt.grid()
                
                # Axes Limit
                ax.set_ylim(1e-17,1e-7)
                
                # Day counter
                text = 'Day: ' + str(days[k])
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                plt.text(0.0465, 0.9, text, # very nice box
                         transform=ax.transAxes, bbox = props)
                
                # Save
                plt.savefig('savedfigs/evo/' + thing + str(k)+ '.png')
                k +=1 
                gc.collect()
# φ - Velocity --------------------------------------------------------------------------            
    if quantity == 'vel':
        thing = coord+quantity
        k = 0
        # fucky_units = True # I forgot to convert to Km/s from m/s
        for part in range(4):
            part = '-' + str(part+1)
            ta_rden = np.load('tavg-data2/evo'+thing+stream[0]+part+ '.npy')
            ta_rden_stream = np.load('tavg-data2/evo'+thing+stream[1]+part+'.npy')
            for i in range(len(ta_rden)):
                # Plot
                fig, ax = plt.subplots(1, clear=True)
                # if fucky_units:
                #     ta_rden[i] = np.multiply(ta_rden[i], 1e-3)
                #     ta_rden_stream[i] = np.multiply(ta_rden_stream[i], 1e-3)
                    
                ax.plot(phis, ta_rden[i], 
                         color = 'maroon', label='All outflows')
                ax.plot(phis, ta_rden_stream[i], 
                         color = 'forestgreen', linestyle='--', label='Excluding the stream')
                ax.set_title(r'Mass-Weighted $\textbf{Outflow}$ velocity for days 21+', fontsize = 18)
                ax.set_xlabel(r'Longitudinal Angle $\left[ \phi \right]$ ', fontsize = 14)
                ax.set_ylabel(r'$|V| \left[ \frac{\mathrm{km}}{\mathrm{s}} \right] $',fontsize = 14)
                ax.set_yscale('log')
                ax.legend(loc='lower left')
                ax.tick_params(axis = 'both', which = 'both', direction='in')
                plt.grid()
                
                # Axes Limit
                ax.set_ylim(5e-2,5e4)
                
                # Light speed
                #ax.axhline(300_000 , c='k')
                # plt.text(0.1, 0.63, 'V = c', # very nice box
                #          transform=ax.transAxes, bbox = props)
                
                # Day counter
                text = 'Day: ' + str(days[k])
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
                plt.text(0.89, 0.1, text, # very nice box
                          transform=ax.transAxes, bbox = props)
                
                # Save
                plt.savefig('savedfigs/evo/' + thing + str(k)+ '.png')
                k +=1 