#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:32:03 2023

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
from astropy.coordinates import cartesian_to_spherical
R, THETA, PHI = cartesian_to_spherical(positions[:,0],
                                       positions[:,1], 
                                       positions[:,2])
# Store it in an array of R, θ, φ, log(ρ)
R = R.value 
THETA = THETA.value
PHI = (-1) * (PHI.value  - np.pi) # (0,2π) -> (-π,π)

# Define new grid
radii = np.linspace(50, 60_000, num = 1_000)
thetas = np.linspace(-np.pi/2, np.pi/2, num = 45) #180)
phis = np.linspace(-np.pi, np.pi, num = 90)

Mass = np.load(fix + '/Mass_' + fix + '.npy')

#%% EVOKE!

mass_grid = THE_TRIPLE_CASTER(radii, R, thetas, THETA, phis, PHI, 
                              Mass)
energy_grid = THE_TRIPLE_CASTER(radii, R, thetas, THETA, phis, PHI, 
                              Orbital)
#%% Isolate unbound
# Do them seperately because we need log scaling and it does not work
# with our current notation.
unbound_mask = np.where(energy_grid > 0, 1, np.nan) # returns 1 on unbound gas, -1
                                                # on bound
unbound_mass = np.multiply(mass_grid, unbound_mask)
unbound_mass = np.log10(unbound_mass)
# unbound_mass = np.nan_to_num(unbound_mass,nan=0, neginf = 0)
custom = plt.cm.jet
custom.set_bad('black',alpha = 1)
# Bound
# bound_mask = np.where(energy_grid < 0, 1, 0) 
# bound_mass = np.multiply(mass_grid, bound_mask)
# bound_mass = np.log10(bound_mass)
# bound_mass = np.nan_to_num(bound_mass, neginf = 0)
# # Make the bounds negative for plotting purposes
# # Recombine
# plot_mass = np.add(bound_mass, unbound_mass) 
#%% Make Plots 
# Can do up to 4300 before restarting kernel
plt.ioff()
import gc
for i in range(len(radii)-1):
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='mollweide')
    # ax2 = fig.add_subplot(122, projection='mollweide')
    
    img1  = ax1.pcolormesh(phis, thetas, unbound_mask[i], cmap = custom,
                           vmin = -11, vmax = -3.2)
    
    fig.suptitle('Unbound Material, 36 days after distruption, R: %.1f' % radii[i])
    # plt.grid(color = 'xkcd:baby shit brown')
    # Colorbars
    plt.colorbar(img1)
    plt.tick_params(
        axis='both',         
        which='both',      
        bottom=False,      
        top=False,         
        labelbottom=False,
        labelleft=False)
    plt.savefig(f'savedfigs/mass_test2{i:03d}.png', dpi=150, format = 'png')
    plt.close('all')
    del fig, ax1, img1 # Needed because otherwise Spyder doesn't release the memory
    gc.collect()
#%% Make movie

command = "ffmpeg -framerate 18 -i 'savedfigs/mass_test2%03d.png' \
         -c:v libx264 -pix_fmt yuv420p -loglevel warning savedmovs/bound_test33.mp4"
os.system(command)
#%% HOW A MOLLWEIDE WORKS IN MATPLOTLIB
fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')
thetas = np.linspace(-np.pi/2, np.pi/2, num = 1_000)
phis = np.linspace(-np.pi, np.pi, num = 1_000)
dummy = np.random.randint(0,10, size=(1000,1000))
ax.pcolormesh(phis, thetas, dummy )