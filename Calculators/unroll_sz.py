#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 10:30:26 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
import numba
from calculators.casters import THE_CASTER
from calculators.stream_identificator import THE_MAX_HOLDER
from calculators.eccentricity import e_calc, ta_calc
#%% Constants
# Need to convert density Msol/Rsol^2 to g/cm
Msol_to_g = 1.989e33
Rsol_to_cm = 6.957e10
converter = Msol_to_g / Rsol_to_cm**2

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2
#%%
# Data Load
fix = '820'
    
X = np.load(fix + '/CMx_' + fix + '.npy')
Y = np.load(fix + '/CMy_' + fix + '.npy')
Z = np.load(fix + '/CMz_' + fix + '.npy')
R = np.sqrt(X**2 + Y**2 + Z**2)
Vx = np.load(fix + '/Vx_' + fix + '.npy')
Vy = np.load(fix + '/Vy_' + fix + '.npy')
Vz = np.load(fix + '/Vz_' + fix + '.npy')
V = np.sqrt(Vx**2 + Vy**2 + Vz**2)

Den = np.load(fix + '/Den_' + fix + '.npy')
# Mass = np.load(fix + '/Mass_' + fix + '.npy')

# Calc true anomalies
positions = np.array((X,Y,Z)).T
velocities = np.array((Vx, Vy, Vz)).T
del X, Y, Vx, Vy, Vz # Memory management

ecc, _ = e_calc(positions, velocities)
true_anomaly = ta_calc(ecc, positions, velocities)
true_anomaly = np.nan_to_num(true_anomaly)
del positions, velocities

# Only care about bound
Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
bound_mask = np.where(Orbital < 0, 1, 0) # returns 1 on bound gas, 0 on unbound
Den = np.multiply(Den,  bound_mask)

# Start at 50->1.7 cause PW pot. End at 30_000 -> 4.5 cause it's reasonable
radii = np.linspace(50, 1000, num = 360)
thetas = np.linspace(0, 2*np.pi, num = 360)


#%% rows are r, columns are θ
density = THE_CASTER(radii, R, thetas, true_anomaly, Den)
density = np.nan_to_num(density) # always useful, argmax no work with NaN
max_holder = THE_MAX_HOLDER(radii, thetas, density)

#%% Plot the
from numpy.polynomial.polynomial import Polynomial as poly
# create the relevant max radii
max_r = []
max_t = []
for i in range(len(max_holder)):
    if radii[int(max_holder[i,1])] < 950:
        max_r.append( radii[int(max_holder[i,1])])
        max_t.append( max_holder[i,2])
        
fig, (ax1, ax2) = plt.subplots(1,2, subplot_kw={"projection": "polar"})
ax1.scatter( max_t, max_r, marker='x', c='k', label='CM of the bound stream')
ax2.scatter( max_t, max_r, marker='x', c='k', label='CM of the bound stream')
fit = poly.fit(max_t, max_r, deg=10)
edge_1 = 2.46 # i dont care close to π
edge_2 = 3.79
relevant_t = np.concatenate( (np.linspace(0,edge_1,num=1_000), np.linspace(edge_2, 2*np.pi,num=1_000)))
ax1.plot(thetas, fit(thetas), c = 'aqua', label = 'Fitted CM orbit')
ax2.scatter(relevant_t, fit(relevant_t), c = 'gold', s=8, marker='h', zorder=3, label = 'Arc Length')
fig.suptitle('Parametrizing the CM orbit')
ax1.legend(loc = 'lower center')
ax2.legend(loc = 'lower center')

#%% Unroller
from scipy.integrate import romberg
full_arclength = romberg(fit, 0,edge_1 ) + romberg(fit, edge_2, 2*np.pi) 

fig, (ax1, ax2) = plt.subplots(1,2, subplot_kw={"projection": "polar"})
ax1.scatter( max_t, max_r, marker='x', c='k', label='CM of the bound stream')
ax1.scatter(relevant_t, fit(relevant_t), c = 'gold', s=8, marker='h', zorder=3, label = 'Arc Length')
ax2.scatter( max_t, max_r, marker='x', c='k', label='CM of the bound stream')
ax2.scatter(relevant_t, fit(relevant_t), c = 'gold', s=8, marker='h', zorder=3, label = 'Arc Length')

test1_arclength = romberg(fit, 1, edge_1)
test1_range = np.linspace(1,edge_1,num=1_000)
ax1.scatter(test1_range, fit(test1_range), c = 'purple', s=8, marker='h', zorder=3)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax1.text(1.85, 700, 's = ' + str(np.round(test1_arclength/full_arclength,3)), bbox=props)

test2_arclength = romberg(fit, 0, edge_1) + romberg(fit, 4, 2*np.pi)
test2_range = np.concatenate((np.linspace(0,edge_1,num=1_000), np.linspace(4,2*np.pi,num=1_000))) 
ax2.scatter(test2_range, fit(test2_range), c = 'teal', s=8, marker='h', zorder=3)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax2.text(1.85, 700, 's = ' + str(np.round(test2_arclength/full_arclength,3)), bbox=props)
fig.suptitle('Examples of parametrization')

def arc(fit,thetas):
    edge_1 = 2.46 # i dont care close to π
    edge_2 = 3.79
    arcs = np.zeros(len(thetas))
    current_progress = 0
    for i in range(len(thetas)):
        if thetas[i] < edge_1:
            arcs[i] = romberg(fit, thetas[i], edge_1)
        elif thetas[i] > edge_2:
            arcs[i] = romberg(fit,0, edge_1) + romberg(fit, thetas[i], 2*np.pi)
        else:
            arcs[i] = romberg(fit, 0, edge_1)
        
        # Progress check
        progress = int(np.round(i/len(thetas),1) * 100)
        if i % 100  == 0 and progress != current_progress:
            print('Arc Calculations', progress, '% done')
            current_progress = progress
    return arcs
#%% Now we cast
fullarc = romberg(fit, 0,edge_1 ) + romberg(fit, edge_2, 2*np.pi) 
inv_fullarc = 1/fullarc
# We have Z
zs = np.linspace(-250, 250, num=100)
# Calc S
S = arc(fit,true_anomaly)
S = np.multiply(S, inv_fullarc)
ss = np.linspace(0,1, num=100)
Entropy = np.load(fix + '/Entropy_' + fix + '.npy')
E_cast = THE_CASTER(zs, Z, ss, S, Entropy)