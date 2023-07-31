#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 16:35:11 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [8 , 4]
import numba
import matplotlib.colors as colors
import colorcet as cc # cooler colormaps
from calculators.casters import THE_CASTER
from scipy import ndimage # rotates image
from extractors.time_extractor import days_since_distruption
from calculators.entropy_stream_identificator import entropy_id_vec_with_neg

# Choose snapshot
fix = '820'
days = str(np.round(days_since_distruption(fix+'/snap_'+fix+'.h5'),1))
# CM Position Data
X = np.load(fix + '/CMx_' + fix + '.npy')
Y = np.load(fix + '/CMy_' + fix + '.npy')

# Import Density
Den = np.load(fix + '/Den_' + fix + '.npy')

# Import Entropy, Mask the stream
Entropy = np.load(fix + '/Entropy_' + fix + '.npy')
stream_mask = entropy_id_vec_with_neg(Entropy) # 1 is stream, -1 is no stream
#%%
# Need to convert Msol/Rsol^2 to g/cm
Msol_to_g = 1.989e33
Rsol_to_cm = 6.957e10
converter = Msol_to_g / Rsol_to_cm**2
Den *=  converter

# Specify new grid:
x_start = -15_000
x_stop = 2000
x_num = 300 # np.abs(x_start - x_stop)
xs = np.linspace(x_start, x_stop, num = x_num )
y_start = -4000
y_stop = 4000
y_num = 300 # np.abs(y_start - y_stop)
ys = np.linspace(y_start, y_stop, num = y_num)

# EVOKE
den_cast = THE_CASTER(xs, X, ys, Y, Den)
stream_mask_cast = THE_CASTER(xs, X, ys, Y, stream_mask)

#%% Remove bullshit and fix things
den_cast = np.nan_to_num(den_cast.T)
den_cast = np.log10(den_cast) # we want a log plot
den_cast = np.nan_to_num(den_cast, neginf=0) # fix the fuckery

# stream_mask_cast = np.nan_to_num(stream_mask_cast.T)
stream_mask_cast[stream_mask_cast>0] = 1
stream_mask_cast[stream_mask_cast<0] = -1
den_cast = np.multiply(den_cast, stream_mask_cast.T)
#%%
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
        
    # den_cast = ndimage.rotate(den_cast, 0)

img = plt.pcolormesh(xs, ys, pixels, cmap='cet_bky')
plt.colorbar(img)
plt.title(r'XY density ($\log_{10}( \rho )$) projection for : '+ days + ' days after disruption')
plt.xlabel('X-coordinate $[R_{\odot}$]')
plt.ylabel('Y-coordinate $[R_{\odot}$]')



