#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 17:02:16 2023

@author: konstantinos
"""

# Vanilla imports
import numpy as np
import gc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# Pretty plots
import matplotlib as mpl
import matplotlib.patheffects as PathEffects
mpl.rcParams.update(mpl.rcParamsDefault)

# plt.rcParams['text.usetex'] = False
import colorcet as cc # cooler colormaps
# Custom Imports
from casters import THE_CASTER
from scipy import ndimage # rotates image
from time_extractor import days_since_distruption
from entropy_stream_identificator import entropy_id_vec_with_neg

fixes = np.arange(750, 1008+1)
folder = 'tde_data/'
for fix in fixes:
    fix = str(fix)
    # Choose snapshot
    days = str(np.round(days_since_distruption('tde_data/snap_'+fix+'/snap_'+fix+'.h5'),1))
    
    # CM Position Data
    X = np.load(folder + 'snap_'+fix + '/CMx_' + fix + '.npy')
    Y = np.load(folder + 'snap_'+fix + '/CMy_' + fix + '.npy')
    
    # Import Density
    Den = np.load(folder + 'snap_'+fix + '/Den_' + fix + '.npy')
    
    # Import Entropy, Mask the stream
    Entropy = np.load(folder + 'snap_'+fix+'/Entropy_' + fix + '.npy')
    stream_mask = entropy_id_vec_with_neg(Entropy) # 1 is stream, -1 is no stream
    
    # Need to convert Msol/Rsol^2 to g/cm
    Msol_to_g = 1.989e33
    Rsol_to_cm = 6.957e10
    
    converter = Msol_to_g / Rsol_to_cm**2
    Den *=  converter
    
    # Specify new grid:
    x_start = -15_000
    x_stop = 200
    x_num = 300 # np.abs(x_start - x_stop)
    xs = np.linspace(x_start, x_stop, num = x_num )
    y_start = -2000
    y_stop = 2000
    y_num = 300 # np.abs(y_start - y_stop)
    ys = np.linspace(y_start, y_stop, num = y_num)
    
    	# EVOKE
    den_cast = THE_CASTER(xs, X, ys, Y, Den)
    stream_mask_cast = THE_CASTER(xs, X, ys, Y, stream_mask)
    
    	# Remove bullshit and fix things
    den_cast = np.nan_to_num(den_cast.T) # need the T here
    den_cast = np.log10(den_cast) # we want a log plot
    den_cast = np.nan_to_num(den_cast, neginf=0) # fix the fuckery
    
    # Renormalize the mask and use it
    stream_mask_cast[stream_mask_cast>0] = 1
    stream_mask_cast[stream_mask_cast<0] = -1
    den_cast = np.multiply(den_cast, stream_mask_cast.T) # need the T here
    
    # Color re-normalization
    pixels = np.zeros_like(den_cast)
    for i in range(len(den_cast)):
        for j in range(len(den_cast)):
            pixel = den_cast[i,j]
            pixels[i,j] = pixel
        
        if pixel > 0 and pixel < 0.1:
            pixels[i,j] = 0
        if pixel > 6:
            pixels[i,j] = 6
        if pixel < 0 and pixel > -0.1:
            pixels[i,j] = 0
        if pixel < -6:
            pixels[i,j] = -6
    
    	# Plotting
    mpl.rcParams.update(mpl.rcParamsDefault)	
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [8 , 4]
    plt.rcParams['font.family'] = 'serif'
    fig = plt.figure(num=1, clear=True)
    ax = fig.add_subplot()
    img = ax.pcolormesh(xs, ys, pixels, cmap='cet_bky')
    	
    	
    fig.colorbar(img)
    ax.set_title('XY Density log10(rho) [g/cm2] projection | Positive is stream') 
    txt = ax.text(0.1, 0.1, 'Day: ' + days,
    color='white', fontweight = 'bold', fontfamily = 'monospace', transform=ax.transAxes)
    txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='k')])
    ax.set_xlabel('X-coordinate [Rsol]')
    ax.set_ylabel('Y-coordinate [Rsol]') 
    plt.savefig('figs/6/stream-check-2/'+fix+'.png')
    	# plt.close(fig)
    gc.collect()
