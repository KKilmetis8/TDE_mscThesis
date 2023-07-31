#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 10:50:45 2023

@author: konstantinos
dai3 for alice FULL
"""
import numpy as np
import numba
from casters import THE_CASTER
from time_extractor import days_since_distruption
from astropy.coordinates import cartesian_to_spherical
from entropy_stream_identificator import entropy_id_vec

m = 4 # 4 or 6

@numba.njit
def time_averager(evolution):
    time_averaged = np.zeros( np.shape(evolution[0]) )
    for i in range(len(evolution)): # Loop over time
        time_averaged = np.add(time_averaged, evolution[i])
    time_averaged = np.multiply(time_averaged, 1/len(evolution))
    return time_averaged

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 10**m # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2

Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
apocenter = 2 * Rt * Mbh**(1/3)

# Need to convert velocity Rsol/time to km/s
Rsol_to_km = 6.957e5
converter_vel = Rsol_to_km / t 

# Need to convert density Msol/Rsol^3 to g/cm^3
Msol_to_g = 1.989e33
Rsol_to_cm = 6.957e10
converter_den = Msol_to_g / Rsol_to_cm**3


#%%
# Batch fixes so it's one call.
if m == 6:
    part = 0
    fixes1 = np.arange(820, 900)
    fixes2 = np.arange(900, 970)
    fixes3 = np.arange(970,1008)
    fixes_list = [fixes1, fixes2, fixes3]
    do_i_care_about_the_stream = [False, True]
    pre = 'tde_data/snap_'

if m == 4:
    part = 0
    fixes1 = np.arange(220, 300)
    fixes2 = np.arange(301, 380)
    fixes3 = np.arange(381, 412+1)
    fixes_list = [fixes1, fixes2, fixes3]
    do_i_care_about_the_stream = [False]
    pre = 'tde_data2/snap_'
    
for fixes in fixes_list:
    part += 1

    dai_3_den = []
    dai_3_den_stream = []
    dai_3_vel = []
    dai_3_vel_stream = []
    for stream in do_i_care_about_the_stream:
        for fix in fixes:
            # Progress Check
            fix = str(fix)
            day = str(np.round(days_since_distruption(pre+fix+'/snap_'+fix+'.h5'),1))
            print('Day:', day)
            # CM Position Data
            X = np.load(pre + fix + '/CMx_' + fix + '.npy')
            Y = np.load(pre + fix + '/CMy_' + fix + '.npy')
            Z = np.load(pre + fix + '/CMz_' + fix + '.npy')
            R = np.sqrt(X**2 + Y**2 + Z**2)
            # Velocity Data
            Vx = np.load(pre + fix + '/Vx_' + fix + '.npy')
            Vy = np.load(pre + fix + '/Vy_' + fix + '.npy')
            Vz = np.load(pre + fix + '/Vz_' + fix + '.npy')
            V = np.sqrt(Vx**2 + Vy**2 + Vz**2)
            del Vx, Vy, Vz
            # Density Data
            Den = np.load(pre + fix + '/Den_' + fix + '.npy')
            Mass = np.load(pre + fix + '/Mass__' + fix + '.npy')
            R, THETA, PHI = cartesian_to_spherical(X,Y,Z)
            R = R.value
            THETA = THETA.value
            
            if stream:
                # Import Entropy, Mask the stream
                Entropy = np.load(pre + fix + '/Entropy_' + fix + '.npy')
                stream_mask = entropy_id_vec(Entropy) # Identify the stream
                stream_mask = (stream_mask - 1) * (-1) # 1 where there is NO stream
                
            # We only care about unbound material
            Orbital = (0.5 * V**2 ) - 10**m / (R-rg)
            unbound_mask = np.where(Orbital > 0, 1, 0) # returns 1 on unbound gas, 0 on bound
            
            # Convert to km/s, g/cm3 and mass weigh
            V *= converter_vel 
            Den *= converter_den
            
            # Combine the masks
            if stream:
                mask = np.multiply(stream_mask, unbound_mask)
                Unbound_V = np.multiply(V, mask)
                Unbound_Den = np.multiply(Den, mask)
            else:
                Unbound_V = np.multiply(V, unbound_mask)
                Unbound_Den = np.multiply(Den, unbound_mask)

            # Specify new grid:
            r_start = 0.2 * 2 * Rt
            r_stop = apocenter
            r_num = 100 # np.abs(x_start - x_stop)
            radii = np.linspace(r_start, r_stop, num = r_num)
            theta_num = 7
            thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
            
            # EVOKE
            v_cast = THE_CASTER(radii, R, thetas, THETA, Unbound_V, 
                                weights = Mass, avg=True)
            v_cast = np.nan_to_num(v_cast)
            d_cast = THE_CASTER(radii, R, thetas, THETA, Unbound_Den, 
                                weights = Mass, avg=False)
            d_cast = np.nan_to_num(d_cast)
            
            if stream:
                dai_3_vel_stream.append(v_cast)
                dai_3_den_stream.append(d_cast)
            else:
                dai_3_vel.append(v_cast)
                dai_3_den.append(d_cast)
    
    np.save('tavg4/' + str(m) + '-bins-rho'+str(part), dai_3_den)
    np.save('tavg4/' + str(m) + '-bins-vel'+str(part), dai_3_vel)
    np.save('tavg4/' + str(m) + '-bins-rho-stream'+str(part), dai_3_den_stream)
    np.save('tavg4/' + str(m) + '-bins-vel-stream'+str(part), dai_3_vel_stream)
