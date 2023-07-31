#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:35:38 2023

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
V = np.sqrt(Vx**2 + Vy**2 + Vz**2) * (Rsol/t) * 1e-3 # km/s
velocities = np.array((Vx, Vy, Vz)).T

# Calculate orbital energy Kin+PW Pot
Orbital = (0.5 * V**2 ) - Mbh / (R-rg)

from astropy.coordinates import cartesian_to_spherical
R, THETA, PHI = cartesian_to_spherical(positions[:,0],
                                       positions[:,1], 
                                       positions[:,2])
# Store it in an array of R, θ, φ, log(ρ)
R = R.value 
THETA = THETA.value
PHI = PHI.value  - np.pi

# Define new grid
radii = np.logspace(1.5, 4.7,num = 1_000)
thetas = np.linspace(-np.pi/2, np.pi/2, num = 180)
phis = np.linspace(-np.pi, np.pi, num = 360)

#%% EVOKE!
vel_grid = THE_TRIPLE_CASTER(radii, R, thetas, THETA, phis, PHI, 
                              V)
vel_grid = np.nan_to_num(vel_grid)

energy_grid = THE_TRIPLE_CASTER(radii, R, thetas, THETA, phis, PHI, 
                              Orbital)

unbound_mask = np.where(energy_grid > 0, 1,0) # returns 1 on unbound gas
# Apply masks
vel_grid = np.multiply(vel_grid, unbound_mask)

#%% Make Plots
plt.ioff()
import gc
for i in range(len(radii)):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='mollweide')
    img = ax.pcolormesh(phis, thetas, vel_grid[i], cmap = 'cet_rainbow4')
    ax.set_title(r'$V$ [km/s], 36 days after \
                 distruption, R: %.1f' % radii[i])
    # plt.grid(color = 'xkcd:baby shit brown')
    plt.colorbar(img)
    plt.tick_params(
        axis='both',         
        which='both',      
        bottom=False,      
        top=False,         
        labelbottom=False,
        labelleft=False)
    plt.savefig(f'savedfigs/Density_{i:03d}.png', dpi=150, format = 'png')
    plt.close('all')
    del fig, ax # Needed because otherwise Spyder doesn't release the memory
    gc.collect()
#%% Make movie

command = "ffmpeg -framerate 28 -i 'savedfigs/Density_%03d.png' \
         -c:v libx264 -pix_fmt yuv420p -loglevel warning savedmovs/den_mov.mp4"
os.system(command)