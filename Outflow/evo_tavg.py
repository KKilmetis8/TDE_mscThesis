#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 19:04:18 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
# Pretty Plots
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8 , 4]
plt.rcParams['axes.facecolor']='whitesmoke'
AEK = '#F1C410' # Important color

# Specs
coord = 'r-'
quantity = 'vel' # vel den
name4 = 'products/evo-M4' + coord + quantity + 'no-'
name6 = 'products/evo-M6' + coord + quantity + 'no-'
# Constants 6
Mbh4 = 1e4
t_fall4 = 40 * (Mbh4/1e6)**(0.5) # days EMR+20 p13
Rt4 =  Mbh4**(1/3) # Msol = 1, Rsol = 1
apocenter4 = 2 * Rt4 * Mbh4**(1/3)
parts4 = 2
tfb4 = np.load('tfb4.npy')
start4 = np.where(tfb4==0.5)
start4 = start4[0][0]
mid4 = np.where(tfb4==0.9975)
mid4 = mid4[0][0] - start4
# 4
Mbh6 = 1e6
t_fall6 = 40 * (Mbh6/1e6)**(0.5) # days EMR+20 p13
Rt6 =  Mbh6**(1/3) # Msol = 1, Rsol = 1
apocenter6 = 2 * Rt6 * Mbh6**(1/3)
parts6 = 4
tfb6 = np.load('tfb6.npy')
start6 = np.where(tfb6==0.5)
start6 = start6[0][0]
mid6 = np.where(tfb6==1.00325)
mid6 = mid6[0][0] - start6


# Data load
for i in range(parts4):
    part = str(i + 1)
    temp4 = np.load(name4 + part + '.npy')
    if i==0:
        data4 = temp4
    else:
        data4 = np.concatenate((data4, temp4))
for i in range(parts6):
    part = str(i + 1)
    temp6 = np.load(name6 + part + '.npy')
    if i==0:
        data6 = temp6
    else:
        data6 = np.concatenate((data6, temp6))        
            
#%%
def time_averager(evolution):
    shape = np.shape(evolution)
    time_averaged = np.zeros(shape[1])
    for i in range(shape[0]): # Loop over time
        time_averaged = np.add(time_averaged, evolution[i])
        
    time_averaged = np.multiply(time_averaged, 1/shape[0])
    return time_averaged

# Tavg 4
data_early4 = data4[:mid4]
data_late4 = data4[mid4:]
data_early4 = time_averager(data_early4)
data_late4 = time_averager(data_late4)
# Tavg 6
data_early6 = data6[:mid6]
data_late6 = data6[mid6:]
data_early6 = time_averager(data_early6)
data_late6 = time_averager(data_late6)
#%%
fig, ax1 = plt.subplots()

    
ax2 = ax1.twinx()


log = True
if log:
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    
if coord == 'r-':
    # New grid
    r_start = 0.5*Rt4
    r_stop = apocenter4
    r_num = 200
    new_grid = np.linspace(r_start, r_stop, num = r_num) / apocenter4   
    # Label
    ax1.set_xlabel('Radial Coordinate $[r/R_a]$')
#
if coord == 'theta-':
    theta_num = 200
    new_grid = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
   
    # Label
    ax1.set_xlabel(r'Latitudinal Angle $\left[ \theta \right]$')
#
if coord == 'phi-':
    phi_num = 200
    new_grid = np.linspace(-np.pi, np.pi, num = phi_num)
    
    ax1.set_xlabel(r'Longitudinal Angle $\left[ \phi \right]$')
#
if quantity == 'den':
    # Label
    ax1.set_ylabel('$10^6 M_\odot$ Density $[g/cm^3]$')
    ax2.set_ylabel(r'$10^4 M_\odot$ Density $[g/cm^3]$', 
                   rotation = 270, labelpad=15)
#    
if quantity == 'vel':
    # Label
    ax1.set_ylabel(r'$10^6 M_\odot$ Velocity $[km/s]$')
    ax2.set_ylabel(r'$10^4 M_\odot$ Velocity $[km/s]$', 
                   rotation = 270, labelpad=15)


ax1.plot(new_grid, data_early6, color = 'k', label = 'Early $10^6 M_\odot$')
ax1.plot(new_grid, data_late6, color = 'k', label = 'Late $10^6 M_\odot$', linestyle = 'dashed')
ax2.plot(new_grid, data_early4, color = AEK, label = 'Early $10^4 M_\odot$')
ax2.plot(new_grid, data_late4, color = AEK, label = 'Late $10^4 M_\odot$', linestyle = 'dashed')

# Legend 
ax1.plot([], [], color = AEK, label = 'Early $10^4 M_\odot$')
ax1.plot([], [], color = AEK, label = 'Late $10^4 M_\odot$', linestyle = 'dashed')
ax1.legend()
