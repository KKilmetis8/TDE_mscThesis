#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 11:28:20 2023

@author: konstantinos

Ypsilon test
"""

import numpy as np

# Custom Imports
from src.Calculators.reducers import THE_REDUCER_2, THE_FILTERER, THE_PIXELATOR
from src.Extractors.time_extractor import days_since_distruption

pre1 = 'new'
check1 = 'tde_data2' + pre1 + '/'
fixes1 = np.arange(205, 207)

pre2 = 'hr4'
check2 = 'tde_data2' + pre2 + '/'
fixes2 = np.arange(204, 206)

def gridder(check, fix, project = False):
    m = 4
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    apocenter = 2 * Rt * Mbh**(1/3)
    fix = str(fix)
    name = check + 'snap_' + fix + '/'
    
    days = str(np.round(days_since_distruption(name + '/snap_' + fix + '.h5'),3))
    # CM Position Data
    X = np.load(name +'CMx_' + fix + '.npy')
    Y = np.load(name +'CMy_' + fix + '.npy')
    Z = np.load(name +'CMz_' + fix + '.npy')
    # Import Density
    Den = np.load(name +'Den_' + fix + '.npy')
    Vol = np.load(name +'Volume_' + fix + '.npy')

    # Need to convert Msol/Rsol^2 to g/cm
    Msol_to_g = 1.989e33
    Rsol_to_cm = 6.957e10
    converter = Msol_to_g / Rsol_to_cm**2
    Den *=  converter
    
    # Specify new grid:
    x_start = -400 # -apocenter - 4 *2*Rt
    x_stop = 100 # 10 * 2*Rt
    x_num = 200 # np.abs(x_start - x_stop)
    xs = np.linspace(x_start, x_stop, num = x_num )
    # y +- 150, z +- 50
    y_start = -200 # -apocenter
    y_stop = 200 # apocenter
    y_num = 100 # np.abs(y_start - y_stop)
    # x_start = -apocenter - 4 *2*Rt
    # x_stop = 10 * 2*Rt
    # x_num = 200 # np.abs(x_start - x_stop)
    # xs = np.linspace(x_start, x_stop, num = x_num )
    # # y +- 150, z +- 50
    # y_start = -apocenter
    # y_stop = apocenter
    # y_num = 100 # np.abs(y_start - y_stop)
    
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
    
    if project:
        den_project = np.zeros( (shape[0], shape[1]) )
        for i in range(len(zs)):
            den_project = np.add(den_project, den_cast[...,i])
    
        return den_project, xs, ys, days
    
    return den_cast, xs, ys, days

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
    
def stacker(check, fixes):
    # First
    den, xs, ys, day_start = gridder(check, fixes[0], project = True)
    
    # Day holder
    days = [day_start]
       
    # Rest 
    for i in range(1, len(fixes)):
        den_new, _, _, day = gridder(check, fixes[i], project = True)
        den = np.add(den, den_new)
        days.append(day)

    # Mean
    inv_total_fixes = 1/len(fixes)
    den = np.multiply(den, inv_total_fixes)
    
    return den, xs, ys, days
    
# Do the thing
den_baseline, xs, ys, days1 = stacker(check1, fixes1)
den_check, _, _, days2 = stacker(check2, fixes2)
    
ypsilon = ypsilon_maker(den_baseline, den_check)
den_baseline = plot_prepare(den_baseline)
den_check = plot_prepare(den_check)

time = str(days1[0]) + '-' + str(days1[-1])

# Save to file
np.save('products/convergance/ypsilon-' + pre1 + '-' + pre2 + '-' + time, ypsilon)
np.save('products/convergance/proj-' + pre1 + '-' + time, den_baseline)
np.save('products/convergance/proj-' + pre2 + '-' + time, den_check)