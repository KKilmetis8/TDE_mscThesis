# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 19:17:59 2023

@author: Konstantinos
"""
# Imports and constants
import numpy as np
import matplotlib.pyplot as plt
import numba
from src.Extractors.time_extractor import linear_fit_days
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1

# Boltzmann Constant in Sim Units
kb_SI = 1.38e-23 # [m^2 kg s^-2 K^-1]
kb = kb_SI * t**2 / (Rsol**2 * Msol) # Kelvins are sim units
#%% Functions we need
#@numba.jit(nopython=True)
def mask_array(arr, mask):
    ''' 
    Numba doesn't support masking so this has to be done manually
    '''
    filtered_arr = np.zeros(len(mask))
    k = 0
    for i in mask:
        filtered_arr[k] = arr[i]
        k = k+1
    return filtered_arr

#@numba.jit(nopython=True)
def cold_masker(rate, threshold):
    '''
    Seperates the dynamically cold stream from the rest of the background.

    Parameters
    ----------
    rate : array-like
        Thermal energy over Orbital enregy.
    threshold : float
        Threshold of coldness. Anything colder will be classified as part of
        the stream.

    Returns
    -------
    cold_enough : array
        Array of indices of cold enough debris. Use as a mask.

    '''
    cold_enough = []
    for i in range(len(rate)):
        if rate[i] < threshold:
            cold_enough.append(i)
    return cold_enough

def ef3_loader(fix, stream=False):
    '''
    Loads data into enclosed_mass because enclosed_mass is numba accelerated
    and can't do it by itself. A wrapper, more or less.
    
    Parameters
    ----------
    fix : str,
        The snapshot for which we want to generate the numbers
        for extended figure 3.
    
    stream: bool, default False
          If True it considers only the dynamically cold stream for
          the calculation.

    Returns
    -------
    enclosed_masses : array,
        Mass enclosed in 60, 100, 200, 300, 1000 Rsol, in Msol

    '''
    # Constants
    G = 6.6743e-11 # SI
    Msol = 1.98847e30 # kg
    Rsol = 6.957e8 # m
    t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
    Mbh = 1e6 # * Msol
    # Need these for the PW potential
    c = 3e8 * t/Rsol # c in simulator units.
    rg = 2*Mbh/c**2
    
    X = np.load(fix+'/CMx_'+fix+'.npy')
    Y = np.load(fix+'/CMy_'+fix+'.npy')
    Z = np.load(fix+'/CMz_'+fix+'.npy')
    T = np.load(fix+'/T_'+fix+'.npy')
    Vx = np.load(fix+'/Vx_'+fix+'.npy')
    Vy = np.load(fix+'/Vy_'+fix+'.npy')
    Vz = np.load(fix+'/Vz_'+fix+'.npy')
    mass = np.load(fix+'/Mass_'+fix+'.npy')
    
    R = np.sqrt( np.power(X,2) + np.power(Y,2)+ np.power(Z,2))
    V = np.sqrt( np.power(Vx,2) + np.power(Vy,2)+ np.power(Vz,2))
    kin = 0.5 * V**2 
    pot = - Mbh / (R-rg)
    enclosed_masses = enclosed_mass(X, Y, Z, T, pot, kin, mass, stream)
    return enclosed_masses

#@numba.jit(nopython=True)
def enclosed_mass(X, Y, Z, T, pot, kin, mass, stream):
    '''
    Calculates the enclosed mass in a few different radii
    (60, 100, 200, 300, 1000 ) away from the SMBH. 

    Parameters
    ----------
    X : arr,
        X coordinate of debris.
    Y : arr,
        Y coordinate of debris.
    Z : arr,
        Z coordinate of debris.
    T : arr,
        Temperature of debris.
    pot : arr,
        Potential Energy of Debris.
    kin : arr,
        Kinetic Energy of Debris.
    mass : arr,
        Mass of Debris.
    stream : bool,
        If True it considers only the dynamically cold stream for
        the calculation..

    Returns
    -------
    mass_in_said_rs : arr,
        The enclosed mass.

    '''
    # Calculate r
    r = np.sqrt( X**2 + Y**2 + Z**2)
    
    # Identifies the Stream, if we want to.
    if stream:
        # Calculate energies and their rate
        orbital = np.abs(np.add(pot,kin))
        thermal = np.multiply(T, kb)
        rate = np.divide(thermal,orbital)
        
        # Get the mask
        cold_enough = cold_masker(rate, 6e-62) # Arbitrary threshold
        
        # Use the mask
        r = r[cold_enough] # mask_array(r, cold_enough)
        mass = mass[cold_enough] # mask_array(mass, cold_enough)
        
    # If object is within the cirical radii, keep its mass
    rs_we_care_about = [60, 100, 200, 300, 1000]
    mass_in_said_rs = []
    for criterion_r in rs_we_care_about:
        m = 0
        for i in range(len(r)):
            # Enclosed Mass
            if r[i] < criterion_r:
                m = m + mass[i] 
        mass_in_said_rs.append(m)
        
    return mass_in_said_rs
#%% Doing the Thing
# ΝΟΤΕ: 1008 does not work, RAM explodes
fixes = ['820','881', '925', '950']
fixdays = []
r60 = []
r100 = []
r200 = []
r300 = []
r1000 = []
for fix in fixes:
    em = ef3_loader(fix, stream = False)
    fixdays.append(str(linear_fit_days(int(fix))))
    r60.append(em[0])
    r100.append(em[1])
    r200.append(em[2])
    r300.append(em[3])
    r1000.append(em[4])
#%% Plotting
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
# Plot size
plt.rcParams['figure.figsize'] = [5.0, 5.0]
plt.plot(fixdays, r60, 'o-', c='navy', label = 'r = 60 $R_{\odot}$ ')
plt.plot(fixdays, r100,'o-', c='maroon', label = 'r = 100 $R_{\odot}$ ')
plt.plot(fixdays, r200,'o-', c='black', label = 'r = 200 $R_{\odot}$ ')
plt.plot(fixdays, r300, 'o-',c='purple', label = 'r = 300 $R_{\odot}$ ')
plt.plot(fixdays, r1000,'o-', c='darkturquoise', label = 'r = 1000 $R_{\odot}$ ')
plt.yscale('log')
plt.grid()
plt.legend()
plt.title('Enclosed Mass in Multiple radii', fontsize=15, pad=5)
plt.xlabel('Days since Distruption')
plt.ylabel(r'Mass inside r $\left[ M_{\odot} \right]$')

    

    