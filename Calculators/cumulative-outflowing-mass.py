#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:11:41 2023

@author: konstantinos
"""
# Vanilla imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
import colorcet as cc
import numba
import pandas as pd
# Homebrew imports
from calculators.entropy_stream_identificator import entropy_stream_id
from calculators.eccentricity import e_calc, ta_calc
from calculators.casters import THE_CASTER, THE_TRIPLE_CASTER
import os

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2

#%% Data Load

# Pick snapshot
fix = '820'

# Formulate positions
X = np.load(fix + '/CMx_' + fix + '.npy')
Y = np.load(fix + '/CMy_' + fix + '.npy')
Z = np.load(fix + '/CMz_' + fix + '.npy')
R = np.sqrt(X**2 + Y**2 + Z**2)
positions = np.array((X, Y, Z)).T

# Formulate Velocities
Vx = np.load(fix + '/Vx_' + fix + '.npy')
Vy = np.load(fix + '/Vy_' + fix + '.npy')
Vz = np.load(fix + '/Vz_' + fix + '.npy')
V = np.sqrt(Vx**2 + Vy**2 + Vz**2)
velocities = np.array((Vx, Vy, Vz)).T

# Calculate orbital energy Kin+PW Pot
Orbital = (0.5 * V**2 ) - Mbh / (R-rg)

#%% 2-D grid for Bound Mass
del X, Y, Z, Vx, Vy, Vz # Memory management
Entropy = np.load(fix + '/Entropy' + fix + '.npy')

# Specify new grid
radii = np.logspace(1.5, 4.7,num = 1_000)#np.linspace(50, 50_000, num = 1_000)
thetas = np.linspace(0, 2*np.pi, num = 1_000)

# Calculate Ecc, TA
ecc, _ = e_calc(positions, velocities)
true_anomaly = ta_calc(ecc, positions, velocities)

# Get stream mask, for said grid
_ , entropy_mask = entropy_stream_id(positions, velocities, Entropy,
                             radii, thetas)
# We want to exclude the stream, turn 1->0, 0->1
entropy_mask = np.add(entropy_mask, -1)
entropy_mask = np.multiply(entropy_mask, -1)
del ecc, Entropy


# Import Mass and project into a grid
Mass = np.load(fix + '/Mass_' + fix + '.npy')
mass_grid = THE_CASTER(radii, R, thetas, true_anomaly, Mass)
energy_grid = THE_CASTER(radii, R, thetas, true_anomaly, Orbital)

# Calculate unbound
unbound_mask = np.where(energy_grid > 0, 1,0) # returns 1 on unbound gas

# Apply masks
mass_grid = np.multiply(mass_grid, unbound_mask)
mass_grid = np.multiply(mass_grid, entropy_mask)

# Plot
fig = plt.figure()
axs = fig.add_subplot(111, projection='polar')
axs.pcolormesh(thetas, radii, mass_grid, norm = colors.LogNorm(),
                   cmap = 'cet_blues')
fig = plt.figure()
axs = fig.add_subplot(111, projection='polar')
axs.pcolormesh(thetas, radii, unbound_mask,
                   cmap = 'Reds')
#%% Attach coordinates
mass_frame = pd.DataFrame(mass_grid,
                          index = radii, columns = thetas)

#%% Unbound Mass
threshold_radius = 20_000 # Rsol
unbound_mass = 0
for radius in radii[radii > threshold_radius]:
    print(radius)
    unbound_mass += mass_frame.loc[radius].sum()



