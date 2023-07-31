#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 19:01:39 2023

@author: konstantinos
"""
# Vanilla Imports
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [4.0, 4.0]
import numba
import matplotlib.colors as colors
import colorcet as cc # cooler colormaps
# Homebrew Imports
from calculators.eccentricity import e_calc, ta_calc
from calculators.casters import THE_CASTER

# Conversion constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
kb = 1.38e-23 * t**2 / (Rsol**2 * Msol) # Boltzmann

# Masker
def masker(entropy, up, down):
    boolean_table = np.zeros_like(entropy)
    for i in range(len(entropy)):
        for j in range(len(entropy)):
            if entropy[i,j] < up and entropy[i,j] > down:
                boolean_table[i,j] = 1
    return boolean_table

def entropy_id_vec(entropy, upper = 1.5e-6, lower = 1e-8):
    mask = np.zeros(len(entropy))
    for i in range(len(mask)):
        if entropy[i]<upper and entropy[i]>lower:
            mask[i] = 1
    return mask

def entropy_id_vec_with_neg(entropy, upper = 1.5e-6, lower = 1e-8):
    mask = np.zeros(len(entropy))
    for i in range(len(mask)):
        if entropy[i]<upper and entropy[i]>lower:
            mask[i] = 1
        else:
            mask[i] = -1
    return mask

#%%,
def entropy_stream_id(positions, velocities, Entropy,
                      radii = np.linspace(50, 30_000, num = 1_000),
                      thetas = np.linspace(0, 2*np.pi, num = 1_000),
                      upper = 1.5e-6, lower = 1e-8):
    
    # Calculate radial coordinate
    R = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)
    
    # Calculate Ecc, TA
    ecc, _ = e_calc(positions, velocities)
    true_anomaly = ta_calc(ecc, positions, velocities)
    
    # Define new grid
    radii = np.linspace(50, 30_000, num = 1_000)
    thetas = np.linspace(0, 2*np.pi, num = 1_000)
    
    # Cast entropy to new grid
    s_cast = THE_CASTER(radii, R, thetas, true_anomaly, Entropy)
    s_cast = np.nan_to_num(s_cast)
    
    # Mask the stream
    mask = masker(s_cast, upper, lower)
    s_mask = np.multiply(s_cast, mask)
    
    return s_mask, mask

# def energy_stream_id(positions, velocities, orbital, T):
#     '''Plain does not work'''
#     # Calculate dynamical Temprature
#     dyn_temp = np.multiply(kb,T)
#     dyn_temp = np.divide(np.abs(orbital),dyn_temp)
    
#     # Calculate radial coordinate
#     R = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)
    
#     # Calculate Ecc, TA
#     ecc, _ = e_calc(positions, velocities)
#     true_anomaly = ta_calc(ecc, positions, velocities)
    
#     # Define new grid
#     radii = np.linspace(50, 30_000, num = 1_000)
#     thetas = np.linspace(0, 2*np.pi, num = 1_000)
    
#     # Cast the dyn_temp to the new grid
#     dyn_temp_cast = THE_CASTER(radii, R, thetas, true_anomaly, dyn_temp)
#     dyn_temp_cast = np.nan_to_num(dyn_temp_cast)
#     return dyn_temp_cast


#%% 
if __name__ == '__main__':
    # Choose snapshot
    fix = '820'

    # CM Position Data
    X = np.load(fix + '/CMx_' + fix + '.npy')
    Y = np.load(fix + '/CMy_' + fix + '.npy')
    Z = np.load(fix + '/CMz_' + fix + '.npy')
    R = np.sqrt(X**2 + Y**2 + Z**2) # To radial
     
    # Velocity Data
    Vx = np.load(fix + '/Vx_' + fix + '.npy')
    Vy = np.load(fix + '/Vy_' + fix + '.npy')
    Vz = np.load(fix + '/Vz_' + fix + '.npy')

    # More data load
    Den = np.load(fix + '/Den_' + fix + '.npy')
    Entropy = np.load(fix + '/Entropy' + fix + '.npy')

    # Calc true anomalies
    positions = np.array((X, Y, Z)).T
    velocities = np.array((Vx, Vy, Vz)).T
    del X, Y, Z, Vx, Vy, Vz # Memory management
    
    # -------------------------------------------------
    s_mask, mask  = entropy_stream_id(positions, velocities, Entropy)
    
    # Calculate Ecc, TA | NOTE: Redundant but this code is never called again
    ecc, _ = e_calc(positions, velocities)
    true_anomaly = ta_calc(ecc, positions, velocities)
    
    # Define new grid
    radii = np.linspace(50, 30_000, num = 1_000)
    thetas = np.linspace(0, 2*np.pi, num = 1_000)
    
    # Cast density to new grid
    den_cast = THE_CASTER(radii, R, thetas, true_anomaly, Den)
    den_cast = np.nan_to_num(den_cast)
    
    # Use the entropy mask on density
    den_mask = np.multiply(den_cast, mask)
    
    # Overplot to check | Note: .T is neccecary since data is R-θ
    # but matplotlib does polar plots in θ-R
    fig = plt.figure()
    axs = fig.add_subplot(111, projection='polar')
    axs.pcolormesh(thetas, radii, den_cast,
                       cmap = 'cet_blues', norm=colors.LogNorm())
    axs.set_yticklabels([''])
    axs.set_title('Real')
    #axs.set_xlim(0,np.pi)
    # Masked
    fig = plt.figure()
    axs = fig.add_subplot(111, projection='polar')
    axs.pcolormesh(thetas, radii, den_mask,
                       cmap = 'cet_blues', norm=colors.LogNorm())
    axs.set_yticklabels([''])
    axs.set_title('Masked')
   #axs.set_xlim(0,np.pi)





