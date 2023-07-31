#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:50:20 2023

@author: konstantinos
logρ - logT
ion tavg
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
from matplotlib import ticker
import matplotlib.patheffects as PathEffects
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8, 5]
plt.rcParams['axes.facecolor']= 	'whitesmoke'
import colorcet
import numba

loadpath = 'products/ion'

# Read data
final4 = 412
frames4 = 187
fudge = 40
start4 = final4 - frames4 + fudge
d4 = None
counter = 0
for i in range(1, frames4 - fudge):
    counter += 1
    npz4 = np.load(loadpath + '/4/' + str(start4 + i) + '.npz')
    if d4 is None:
        d4 = npz4['den']
        t4 = npz4['temp']
    else:
        d4_n = npz4['den']
        t4_n = npz4['temp']
        d4 = np.add(d4,d4_n)
        t4 = np.add(t4,t4_n)

# Average
d4 = np.multiply(d4, 1/counter)
t4 = np.multiply(t4, 1/counter)
Mbh = 10**4
Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
apocenter = 2 * Rt * Mbh**(1/3)
radii_start = 0.2*2*Rt
radii_stop = apocenter
r4 = np.linspace(radii_start, radii_stop, 1000) / apocenter
peri4 = npz4['peri']

# ΞΑΝΑ
final6 = 1008
frames6 = 400
start6 = final6 - frames6

d6 = None
counter = 0
for i in range(1, frames6):
    counter += 1
    npz6 = np.load(loadpath + '/6/' + str(start6 + i) + '.npz')
    if d6 is None:
        d6 = npz6['den']
        t6 = npz6['temp']
    else:
        d6_n = npz6['den']
        t6_n = npz6['temp']
        d6 = np.add(d6,d6_n)
        t6 = np.add(t6,t6_n)

d6 = np.multiply(d6, 1/(frames6))
t6 = np.multiply(t6, 1/(frames6))
Mbh = 10**6
Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
apocenter = 2 * Rt * Mbh**(1/3)
radii_start = 0.2*2*Rt
radii_stop = apocenter
r6 = np.linspace(radii_start, radii_stop, 1000) / apocenter
peri6 = npz6['peri']
#%% Plotting

# Fig init
fig, ax = plt.subplots(1,2, num=1, clear=True, tight_layout = True,
                       sharex = True, sharey = True)

# Image making
img1 = ax[0].scatter(d6, t6,
                     c = r6, cmap = 'cet_bmy',
                     s=20, marker='h', zorder = 3)
ax[0].scatter(d6[peri6], t6[peri6],
                     c = 'r', cmap = 'cet_bmy',
                     s=40, marker='x', zorder = 3)
# plt.colorbar(img1)
end = len(d4) # 100
img2 = ax[1].scatter(d4[0:end], t4[0:end], 
                     c = r4[0:end], cmap = 'cet_bmy',
                     s=20, marker='h', zorder = 3)
ax[1].scatter(d4[peri4], t4[peri4],
                     c = 'r', cmap = 'cet_bmy',
                     s=40, marker='x', zorder = 3)
cax = fig.add_axes([1, 0.065, 0.02, 0.86])
fig.colorbar(img1, cax=cax)

inset = True
if inset:
    l, b, h, w = .565, .35, .4, .145
    inset = fig.add_axes([l, b, w, h])
    inset.tick_params(axis = 'both', which = 'both', direction='in')
    instart = 300
    instop = 1_000
    instep = 1
    inset.scatter(d4[instart:instop:instep], t4[instart:instop:instep],
                c = img2.to_rgba(r4[instart:instop:instep]),
                s=8, marker='h', zorder = 5)
    inset.set_xlim(-6.5, -5.2)
    inset.set_ylim(4.2, 4.5)
    inset.xaxis.set_major_locator(plt.MaxNLocator(2))
    inset.yaxis.set_major_locator(plt.MaxNLocator(2))
    inset.grid()
     
# # Days text
dayx = 0.38 
dayy = 0.8
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
txt1 = fig.text(dayx, dayy, 'Temporal Range (0.3 - 1.6) t/tfb:',
            fontsize = 10,
		    color='k', fontfamily = 'monospace', bbox = props)
# #txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
# txt2 = ax[1].text(dayx, dayy, 't/tfb: ' + str(days4), fontsize = 10,
# 		    color='k', fontfamily = 'monospace', bbox = props,
#                 transform=ax[1].transAxes)
# Ionized text
ionx = 0.69
iony = 0.19
txt1 = ax[0].text(ionx, iony, '50 \% Ionized', fontsize = 12,
		    color='k', fontfamily = 'monospace', rotation = 8,
                transform=ax[0].transAxes)
#txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
ionx = 0.15
iony = 0.15
txt1 = ax[1].text(ionx, iony, '50 \% Ionized', fontsize = 12,
		    color='k', fontfamily = 'monospace', rotation = 8,
                transform=ax[1].transAxes)
# Ionized text 2
ionx = 0.69
iony = 0.28
txt1 = ax[0].text(ionx, iony, '95 \% Ionized', fontsize = 12,
		    color='k', fontfamily = 'monospace', rotation = 8,
                transform=ax[0].transAxes)
#txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
ionx = 0.15
iony = 0.23
txt1 = ax[1].text(ionx, iony, '95 \% Ionized', fontsize = 12,
		    color='k', fontfamily = 'monospace', rotation = 7,
                transform=ax[1].transAxes)
#txt2.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])

# Distance text 
ionx = 1.08
iony = 0.35
#txt1.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='k')])
txt1 = fig.text(ionx, iony, 'Distance to BH [r/R$_a$]', fontsize = 12,
		    color='k', fontfamily = 'monospace', rotation = 270)
# 90 % line
c2 = -1.578e5
c3 = 4.01e-9 # where we go from NR to ER
y = 0.95
c4 = y**2 / (1-y)
T_line = np.linspace(3, 12, num=200)
den_line = np.log10(c3 / c4) + 1.5 * T_line + c2/10**T_line * np.log10( 2.718281) # e

ax[0].plot(den_line, T_line,
         c = 'k', linestyle='dashed',
         zorder=5)    
ax[1].plot(den_line, T_line,
         c = 'k', linestyle='dashed',
         zorder=5)
inset.plot(den_line, T_line,
         c = 'k', linestyle='dashed',
         zorder=5)  
ax[0].fill_between(den_line, T_line, alpha = 0.2,
                color='r', zorder=1)
ax[1].fill_between(den_line, T_line, alpha = 0.2,
                color='r', zorder=1)
inset.fill_between(den_line, T_line, alpha = 0.2,
                color='r', zorder=1)

# 50 % line
y = 0.5
c4 = y**2 / (1-y)
T_line = np.linspace(3, 12, num=200)

den_line = np.log10(c3 / c4) + 1.5 * T_line + c2/10**T_line * np.log10( 2.718281) # e

ax[0].plot(den_line, T_line,
         c = 'k', linestyle='dashed',
         zorder=5)    
ax[1].plot(den_line, T_line,
         c = 'k', linestyle='dashed',
         zorder=5)
inset.plot(den_line, T_line,
         c = 'k', linestyle='dashed',
         zorder=5)
ax[0].fill_between(den_line, T_line, 
                color='#d2e7d6', zorder=1)
ax[1].fill_between(den_line, T_line,
                color='#d2e7d6', zorder=1)
inset.fill_between(den_line, T_line,
                color='#d2e7d6', zorder=1)


# Axis labels
fig.text(0.5, -0.01, r'$log_{10}(\rho)$  [g/cm$^3$]', ha='center', fontsize = 14)
fig.text(-0.02, 0.5, r'$log_{10}(T)$  [K]', va='center', rotation='vertical', fontsize = 14)
# ax[0].set_xlabel(r'$log_{10}(\rho)$  [g/cm$^3$]', fontsize = 12)
# ax[0].set_ylabel(r'$log_{10}(T)$  [K]', fontsize = 12)
#ax[1].set_xlabel(r'$log_{10}(\rho)$  [g/cm$^3$]', fontsize = 12)
#ax[1].set_ylabel(r'$log_{10}(T)$  [K]', fontsize = 12)
ax[0].tick_params(axis = 'both', which = 'both', direction='in')
ax[1].tick_params(axis = 'both', which = 'both', direction='in')

# Axis lims
xlow = -13
xhigh = -4.5
ylow = 3
yhigh = 8.2
ax[0].set_xlim(xlow, xhigh)
ax[0].set_ylim(ylow, yhigh)
ax[1].set_xlim(xlow, xhigh)
ax[1].set_ylim(ylow, yhigh)

# Titles
# fig.suptitle('Ionization of debris', fontsize = 17)
ax[0].set_title('$10^6 M_\odot$', fontsize = 15)
ax[1].set_title('$10^4 M_\odot$', fontsize = 15)
ax[0].grid(zorder = 3)
ax[1].grid(zorder = 3)
        