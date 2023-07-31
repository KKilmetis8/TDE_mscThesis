#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 11:39:32 2023

@author: konstantinos
"""
#%% Intro stuff.

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [8.0, 6.3]
import matplotlib.colors as colors
import matplotlib.cm as cm
import numba
from calculators.casters import THE_CASTER

# Conversion constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
#%% SUMMON THE MAX HOLDER
@numba.njit
def THE_MAX_HOLDER(radii, thetas, density):
    '''
    Holds the INDICES of max. density in each radial ray
    
    Parameters
    ----------
    radii : arr,
        Array of radii we want to cast to.
    thetas: arr,
        Array of true anomalies we want to cast to.
    density : arr,
        Casted down densities.

    Returns
    -------
    max_holder : arr 1x3,
        Contains the maximum density for every true anomaly.
        
        Structure: density, 
                  index of radius of maximum density,
                  theta

    '''
    max_holder = np.zeros((len(thetas),3))
    
    for i in range(len(thetas)):
        # Grab all the densities for a given theta
        den_column = density[:,i]
        
        # Find the index of the maximum
        idx = np.argmax(den_column)
        
        # Store in max holder
        max_holder[i,0] = den_column[idx]
        max_holder[i,1] = idx # index of max radius is more handy.
        max_holder[i,2] = thetas[i]
        
    return max_holder

#%% CONSTRUCT THE WIDTH FINDER
@numba.njit
def THE_WIDTH_FINDER(radii, thetas, density, max_holder, threshold=1.05, 
                     before = False):
    '''
    Parameters
    ----------
    radii : arr,
        Array of radii we want to cast to.
    thetas: arr,
        Array of true anomalies we want to cast to.
    density : arr,
        Casted down densities.
        
    max_holder : arr 1x3,
        Contains the maximum density for every true anomaly.
        Structure: density, 
                  index of radius of maximum density,
                  theta

    Returns
    -------
    width_holder : TYPE
        DESCRIPTION.

    '''
    width_holder = np.zeros(len(thetas))
    for i in range(len(thetas)):
        current_den = max_holder[i,0]
        if current_den < 1e-20:
            continue
        center_radius_idx = int(max_holder[i,1]) # int fixes array fuckery
        width_counter = 0
        # Start expanding in radius
        # go in both directions since this is largely symmetrical
        for j in range(1, len(radii)):
            # Find the radii one index away from the center one
            if center_radius_idx == 0:
                radius_away = center_radius_idx + j
                
                # Get their density to the total
                new_den = density[radius_away, i]
                
                # Check if their density is a considerable part of the total amount
                is_considerable = (current_den + new_den) / current_den
                
                if is_considerable > threshold:
                    current_den += new_den
                    width_counter += 1
                else:
                    width_holder[i] = width_counter
                    break
            else:
                radius_away = center_radius_idx + j
                
                # In both directions
                radius_closer = center_radius_idx - j
                
                # Get their density to the total
                new_den = density[radius_away, i] + density[radius_closer, i]
                
                # Check if their density is a considerable part of the total amount
                is_considerable = (current_den + new_den) / current_den
                
                if is_considerable > threshold:
                    current_den += new_den
                    width_counter += 1
                else:
                    width_holder[i] = width_counter
                    break
        
        return width_holder

# Queries it
@numba.njit
def THE_ONE_WHO_MAKES_A_BOOLEAN_TABLE_FOR_THE_STREAM(max_holder, width_holder):
    stream = np.zeros(( len(max_holder), len(max_holder) ))
    
    for i in range(len(max_holder)):
        # Unpacking + clarity
        starting_idx = max_holder[i,1]
        width = width_holder[i]
        if starting_idx==0:
            # Generate the indexes which contain the stream
            upper_idx = starting_idx + width
            # radii of stream, true anomaly.
            stream[0:upper_idx,i] = 1
        else:
            # Generate the indexes which contain the stream
            upper_idx = starting_idx + width
            lower_idx = starting_idx - width
            # radii of stream, true anomaly.
            stream[lower_idx:upper_idx,i] = 1
        
    return stream

@numba.njit 
def mass_in_stream(thetas, widths, stream_mass, threshold):
    good_radii = np.zeros(len(thetas))
    
    # Find the radius within which lies 50% of the mass for every theta
    for i in range(len(thetas)):
        mass_counter = 0
        cell_counter = 0
        tot_mass = np.sum(stream_mass[:,i]) # mass in that theta
        max_idx = np.argmax(stream_mass[:,i]) # center of stream
        
        # Iterate over widths
        for j in range(widths[i]):
            
            mass_counter += stream_mass[max_idx + j,i] +\
                            stream_mass[max_idx - j,i]
            cell_counter += 1
            
            if mass_counter >= threshold * tot_mass:
                # Cells -> RSun, /2 to get from diameter to radius
                cell_to_rsun = radii[i] - radii[i-1]
                good_radii[i] = cell_counter * cell_to_rsun / 2 
                break
    return good_radii
#%% EVOKE THEM !
if __name__ == '__main__':
    # Data Load
    fix = '820'
    
    # Check if we need to take into account the coordinate change
    if int(fix) < 720:
        before = True
    else:
        before = False
        
    X = np.load(fix + '/CMx_' + fix + '.npy')
    Y = np.load(fix + '/CMy_' + fix + '.npy')
    Z = np.load(fix + '/CMz_' + fix + '.npy')
    R = np.sqrt(X**2 + Y**2 + Z**2)
    Vx = np.load(fix + '/Vx_' + fix + '.npy')
    Vy = np.load(fix + '/Vy_' + fix + '.npy')
    Vz = np.load(fix + '/Vz_' + fix + '.npy')
    from calculators.eccentricity import e_calc, ta_calc
    Den = np.load(fix + '/Den_' + fix + '.npy')
    # Mass = np.load(fix + '/Mass_' + fix + '.npy')
    
    # Define the smaller version of the grid 

    
    if before:
        true_anomaly = np.arctan2(Y,X)
        thetas = np.linspace(-np.pi,np.pi, num=36)
        radii = np.linspace(0.01, np.max(R), num=36)
    else:
        # Start at 50->1.7 cause PW pot. End at 30_000 -> 4.5 cause it's reasonable
        radii = np.linspace(50, 30_000, num = 360)
        thetas = np.linspace(0, 2*np.pi, num=360)
        # Calc true anomalies
        positions = np.array((X,Y,Z)).T
        velocities = np.array((Vx, Vy, Vz)).T
        del X, Y, Vx, Vy # Memory management
        ecc, _ = e_calc(positions, velocities)
        true_anomaly = ta_calc(ecc, positions, velocities)
        true_anomaly = np.nan_to_num(true_anomaly)
        # rows are r, columns are θ
        
    density = THE_CASTER(radii, R, thetas, true_anomaly, Den)
    density = np.nan_to_num(density) # always useful, argmax no work with NaN
    max_holder = THE_MAX_HOLDER(radii, thetas, density)
    
    #%% Do the thing
    threshold = 1.03
    width_holder = THE_WIDTH_FINDER(radii, thetas, density, max_holder,
                                   threshold)
    
    # Width holder has cell number held, convert to radii
    stream = THE_ONE_WHO_MAKES_A_BOOLEAN_TABLE_FOR_THE_STREAM(max_holder,
                                                              width_holder)
    # Mass in stream
    stream_mass = np.multiply(density, stream)
    stream_mass = np.nan_to_num(stream_mass)
    rad50 = mass_in_stream(thetas, width_holder, stream_mass, 
                           threshold = 0.5)
    rad90 = mass_in_stream(thetas, width_holder, stream_mass, 
                           threshold = 0.9)
    
    # Plot fig2.
    plt.figure()
    plt.plot(thetas, rad50, color='maroon', label = '50 \% mass radius')
    plt.plot(thetas, rad90, '--' ,color='teal' , label = '90 \% mass radius')
    # plt.xlim(3.05, 3.4)
    plt.yscale('log')
    plt.grid()
    plt.legend()
    plt.title('Stream radius enclosing different slices of the mass (EF 2)')
    plt.ylabel('Stream Radius [$R_{\odot}$]')
    plt.xlabel('True Anomaly [rad]')
 
    #%% Plot to Check
    # import colorcet as cc
    # fig, axs = plt.subplots(2,4)
    # thresholds = [1.03]
    # # Make titles
    # titles = ['True']
    # for t in thresholds:
    #     titles.append(str(np.round((t-1)*100,4)) + '\%')
    # # Plot it
    # k = 0
    # for i in range(2):
    #     for j in range(4):
    #         axs[i,j].imshow(streams[k], cmap = 'cet_blues')
    #         axs[i,j].set_xlim(15,80)
    #         axs[i,j].set_ylim(0,100)
    #         axs[i,j].set_title(titles[k])
    #         axs[i,j].grid()
    #         k += 1
    # axs[0,0].imshow(density, cmap = 'cet_blues',
    #                 norm = colors.LogNorm())
    # axs[0,0].set_title('Real Stream')
    # fig.suptitle('Stream Identification for different thresholds', fontsize=16)
    
    # #%% another thing
    # plt.rcParams['figure.figsize'] = [4.0, 4.0]
    # for i in range(len(streams)):
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='polar')
    #     ax.pcolormesh(thetas, radii, streams[i].T, cmap='cet_blues')
    #     ax.set_yticklabels([''])
    #     ax.set_xlim(0,np.pi)
    #     ax.set_title(titles[i])






