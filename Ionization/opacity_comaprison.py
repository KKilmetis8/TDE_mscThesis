#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:47:37 2023

@author: konstantinos
"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import colorcet
import numpy as np
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [6.0, 6.0]
plt.rcParams['axes.facecolor']= 	'whitesmoke'
AEK = '#F1C410' # Important color

Msol_to_g = 1.989e33
Rsol_to_cm = 6.957e10
def opacity_sum(rho,T):
    X = 0.7 # x is actually 0.7
    
    constant = 0.64e23
    kt = 0.2 * (1+X)
    k = constant * rho * T**(-7/2) + kt # kramers + thompson
    return k
def opacity_ff(rho,T):
    # cm2 / g
    cm_to_Rsol = 1/Rsol_to_cm
    g_to_Msol = 1/Msol_to_g
    converter = cm_to_Rsol**2 / g_to_Msol 
    constant = 3.6e22 # * converter # cm2 / g -> Rsol2 / Msol

    k = constant * rho * T**(-7/2) * (1 + 0.7)
    return k

def opacity_t(rho,T):
    X = 0.7 # x is actually 0.7

    # Let's just do Thompson
    k = 0.2 * (1+X) 
    return k

temps = np.logspace(3.5,8, num = 1000)
# mean for 1e6 is 1e-9
opac_s = [ opacity_sum(1e-9, temp) for temp in temps]
opac_ff = [ opacity_ff(1e-9, temp) for temp in temps]
opac_t = [ opacity_t(1e-9, temp) for temp in temps]

plt.plot(temps, opac_ff, color=AEK, label = 'Free-Free')
plt.plot(temps, opac_t, c='k', label = 'Thompson')
plt.plot(temps, opac_s, c='maroon', label='Sum')
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.1, 50)
plt.grid()
plt.title('Opacity Comparison')
plt.xlabel('Temperature [K]')
plt.ylabel('Opacity [cm$^2$/g]')
plt.legend()