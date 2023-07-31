#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 18:11:00 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [8.0, 4.0]
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.facecolor']='whitesmoke'
from extractors.time_extractor import days_since_distruption

pre = 'products'
en4 = np.load(pre + '/Energies-4.npy')
en6 = np.load(pre + '/Energies-6.npy')

fig, ax = plt.subplots(1,2, tight_layout=True)

ax[0].plot(en4[0], en4[1], c='tab:purple', label='Orbital')
ax[0].plot(en4[0], en4[2], c='tab:red', label = 'Internal')
ax[0].plot(en4[0], en4[3], c='tab:olive', label = 'Radiation')

ax[1].plot(en6[0], en6[1],  c='tab:purple', label='Orbital')
ax[1].plot(en6[0], en6[2],  c='tab:red', label = 'Internal')
ax[1].plot(en6[0], en6[3],  c='tab:olive', label = 'Radiation')

# lims
ax[0].set_xlim(0.7,1.65)
ax[1].set_xlim(0.7, 1.65)

# Scale
ax[0].set_yscale('log')
ax[1].set_yscale('log')

ax[0].set_title(r'$10^4$ IMBH')
ax[1].set_title(R'$10^6$ SMBH')

# Axis labels
ax[0].set_xlabel(r' Time / Fallback time $\left[ t/t_{FB} \right]$', fontsize = 14)
ax[1].set_xlabel(r' Time / Fallback time $\left[ t/t_{FB} \right]$', fontsize = 14)
ax[0].set_ylabel(r'Total Energy [$M_{\odot} R_{\odot}^2 \widetilde{t}^{-2}$]', fontsize = 14)
ax[0].tick_params(axis = 'both', which = 'both', direction='in')
ax[1].tick_params(axis = 'both', which = 'both', direction='in')

ax[0].grid()
ax[1].grid()
ax[0].legend()
ax[1].legend()


fig.suptitle('Energy evolution', fontsize = 17)