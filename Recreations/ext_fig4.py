# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 16:33:48 2022

@author: Konstantinos
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 09:21:27 2022

@author: Konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from numba import jit
#%% Data Import
fix = '881' # It's both a prefix and a suffix so they cancel out!
X = np.load(fix + '/X_' + fix + '.npy')
Y = np.load(fix + '/Y_' + fix + '.npy')
Z = np.load(fix + '/Z_' + fix + '.npy')
Den = np.load(fix + '/Den_' + fix + '.npy')

#%% Integration Function
# Using Numbas because this takes ages otherwise

@jit(nopython=True)
def integrate(tointegrate, dimension, chunk_size = 2):
    result = np.zeros( (len(dimension),))
    iterating_range = np.arange(len(dimension), step = chunk_size)
    old_done = 0
    if chunk_size == 1:
        for i in iterating_range:
            y = tointegrate[i:i+chunk_size]
            x = dimension[i:i+chunk_size]
            integral = y*x
            result[i:i+chunk_size] = integral
            # Progress Check
            done = int(100*i/iterating_range[-1])
            if done>old_done:
                old_done = done    
                print( done, ' percent done')
        return result
    else:
        for i in iterating_range:
            y = tointegrate[i:i+chunk_size]
            x = dimension[i:i+chunk_size]
            integral = np.trapz(y, x)
            result[i:i+chunk_size] = integral
            # Progress Check
            done = int(100*i/iterating_range[-1])
            if done>old_done:
                old_done = done    
                print( done, ' percent done')
        return result

rhodz = integrate(Den, X, chunk_size=1)
#%% Units
# its in Msol/Rsol^2 (you integrate once!)
# M_sol to g and Rsol to cm
Msol_to_g = 1.989e33
Rsol_to_cm = 6.957e10
rhodz = rhodz * Msol_to_g / Rsol_to_cm**2
# Log for easier scale
rhodz = np.log10(rhodz)
rhodz = np.nan_to_num(rhodz,nan=0, neginf=0)
rhodz[rhodz<0] = 0
#%% colors
colores = np.zeros(len(rhodz))
sizes = np.ones(len(rhodz)) 
colores = rhodz
# colores[colores>5] = 5
# Ensures Black Background
colores[rhodz<0.5] = 0
sizes[rhodz<0.5] = 100_000 
#%% Sorting
import pandas as pd
df = pd.DataFrame({'Y':Y,
                   'Z':Z,
                   'rhodz':rhodz,
                   'color':colores,
                   'size':sizes})

df.sort_values(by='rhodz', ascending=True, inplace=True)
#%% Plotting
import colorcet as cc

fire = cc.cm.fire
jet = cc.cm.rainbow
plt.figure()
stop = len(X)
size = 1_000_000
start = stop-size

plt.scatter(df['Y'][start:stop], df['Z'][start:stop], 
            c=df['color'][start:stop], s=df['size'][start:stop], cmap=fire)
plt.title(r' Projected $\log_{10}(\rho)$ in the YZ plane for Snapshot 881 $[g/cm^2]$')
plt.xticks(rotation = 45)
plt.xlabel('Y - Coordinate [AU]')
plt.ylabel('Z - Coordinate [AU]')
plt.xlim(-500,500)
plt.ylim(-200,200)
plt.colorbar()