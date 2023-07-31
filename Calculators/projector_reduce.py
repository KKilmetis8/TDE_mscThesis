#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 20:45:15 2023

@author: konstantinos
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 11:28:20 2023

@author: konstantinos

Ypsilon test
"""

import numpy as np
import matplotlib.pyplot as plt
# Custom Imports
from src.Calculators.casters_split import THE_CASTER_SPLIT_2
from src.Calculators.reducers import THE_REDUCER_2, THE_FILTERER, THE_PIXELATOR, THE_PROJECTOR
from src.Extractors.time_extractor import days_since_distruption


fixes = ['177', '232', '263']
fixes = ['263']
choice = 'XZ'

def gridder(fix, choice):
    m = 4
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    apocenter = 2 * Rt * Mbh**(1/3)
    fix = str(fix)
    name = fix + '/'

    days = str(np.round(days_since_distruption(name +'snap_' + fix + '.h5'),3))
    # CM Position Data
    X = np.load(name +'CMx_'  + fix +  '.npy')
    Y = np.load(name + 'CMy_' + fix + '.npy')
    Z = np.load(name + 'CMz_' + fix + '.npy')
    # Import Density
    Den = np.load(name + 'Den_' + fix +'.npy')
    Vol = np.load(name + 'Volume_' + fix + '.npy')
    # Mass = np.load(name + 'Mass.npy')
        
    # Need to convert Msol/Rsol^2 to g/cm
    Msol_to_g = 1.989e33
    Rsol_to_cm = 6.957e10
    converter = Msol_to_g / Rsol_to_cm**2
    Den *=  converter
    
    if choice == 'XY':
        # Specify new grid:
        x_start = -apocenter - 4 *2*Rt
        x_stop = 10 * 2*Rt
        x_num = 100 # np.abs(x_start - x_stop)
        xs = np.linspace(x_start, x_stop, num = x_num )
        # y +- 150, z +- 50
        y_start = -apocenter
        y_stop = apocenter
        y_num = 100 # np.abs(y_start - y_stop)
        ys = np.linspace(y_start, y_stop, num = y_num)    
        z_start = -100
        z_stop = 100
        z_num = 10 # np.abs(y_start - y_stop)
        zs = np.linspace(z_start, z_stop, num = z_num)    
        # EVOKE
        pixels = THE_PIXELATOR(xs, ys, zs)
        nX, nY, nZ, nDen = THE_FILTERER(X, Y, Z, Den, Vol)
        shape = (x_num, y_num, z_num)
        den_cast = THE_REDUCER_2(pixels, shape, 
                                 nX, nY, nZ, 
                                 nDen)
        
        # Remove bullshit and fix things
        den_cast = np.nan_to_num(den_cast) # need the T here
        
        den_cast = THE_PROJECTOR(den_cast, zs, choice)
        
        return den_cast, xs, ys, days
    if choice == 'XZ':
        # Specify new grid:
        x_start = -apocenter - 4 *2*Rt
        x_stop = 10 * 2*Rt
        x_num = 300 # np.abs(x_start - x_stop)
        xs = np.linspace(x_start, x_stop, num = x_num )
        # y +- 150, z +- 50
        y_start = -300 # -apocenter
        y_stop = 300 # apocenter
        y_num = 6 # np.abs(y_start - y_stop)
        ys = np.linspace(y_start, y_stop, num = y_num)    
        z_start = -apocenter
        z_stop = apocenter
        z_num = 300 # np.abs(y_start - y_stop)
        zs = np.linspace(z_start, z_stop, num = z_num)    
        # EVOKE
        pixels = THE_PIXELATOR(xs, ys, zs)
        nX, nY, nZ, nDen = THE_FILTERER(X, Y, Z, Den, Vol)
        shape = (x_num, y_num, z_num)
        den_cast = THE_REDUCER_2(pixels, shape, 
                                 nX, nY, nZ, 
                                 nDen)
        
        # Remove bullshit and fix things
        den_cast = np.nan_to_num(den_cast) # need the T here
        
        den_cast = THE_PROJECTOR(den_cast, ys, choice)
        
        return den_cast, xs, zs, days
    if choice == 'YZ':
        # Specify new grid:
        x_start = -apocenter - 4 *2*Rt
        x_stop = 10 * 2*Rt
        x_num = 30 # np.abs(x_start - x_stop)
        xs = np.linspace(x_start, x_stop, num = x_num )
        # y +- 150, z +- 50
        y_start = -apocenter
        y_stop = apocenter
        y_num = 100 # np.abs(y_start - y_stop)
        ys = np.linspace(y_start, y_stop, num = y_num)    
        z_start = -apocenter
        z_stop = apocenter
        z_num = 100 # np.abs(y_start - y_stop)
        zs = np.linspace(z_start, z_stop, num = z_num)    
        # EVOKE
        pixels = THE_PIXELATOR(xs, ys, zs)
        nX, nY, nZ, nDen = THE_FILTERER(X, Y, Z, Den, Vol)
        shape = (x_num, y_num, z_num)
        den_cast = THE_REDUCER_2(pixels, shape, 
                                 nX, nY, nZ, 
                                 nDen)
        
        # Remove bullshit and fix things
        den_cast = np.nan_to_num(den_cast) # need the T here
        
        den_cast = THE_PROJECTOR(den_cast, xs, choice)
        
        return den_cast, ys, zs, days

def plot_prepare(den_cast, color = True):
    # We want a log plot
    den_cast = np.log10(den_cast) 
    
    # Fix the fuckery
    den_cast = np.nan_to_num(den_cast, neginf=0)
    # Color re-normalization
    if color:
        den_cast[den_cast<0.1] = 0
        den_cast[den_cast>8] = 8
    
    # Transpose to look like we're used to
    den_cast = den_cast.T
    return den_cast

def ypsilon_maker(den_baseline, den_check):
    ypsilon = np.divide(den_baseline, den_check)
    ypsilon = plot_prepare(ypsilon, color = False)
    return ypsilon
    
for fix in fixes:
    den_baseline, xs, ys, day1 = gridder(fix, choice)
    den_cast = plot_prepare(den_baseline)
    plt.figure()
    plt.pcolormesh(xs, ys, den_cast)
    np.save(choice + '4-' + fix, den_cast)

