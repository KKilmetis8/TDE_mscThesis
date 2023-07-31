#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 14:46:18 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [4.0, 4.0]
AEK = '#F1C410' # Important color
import os
alice = False

if alice:
    from casters import THE_CASTER
    from time_extractor import days_since_distruption
else:
    from calculators.casters import THE_CASTER
    from extractors.time_extractor import days_since_distruption
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as PathEffects
    plt.rcParams['text.usetex'] = True
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['figure.figsize'] = [9 , 8]
    import colorcet
    

# Constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
 
#%%
def maker(m, fix, alice):
    Mbh = 10**m
    Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
    t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13
    
    c = 3e8 * t/Rsol # c in simulator units.
    rg = 2*Mbh/c**2
    
    # Choose snapshot
    fix = str(fix)

    if alice:
        Mass = np.load(fix + '/Mass__' + fix + '.npy')
        if m==4:
            folder = 'tde_data2/snap_' + fix
        if m==6:
            folder = 'tde_data/snap_' + fix
    else:
        Mass = np.load(fix + '/Mass_' + fix + '.npy')
        folder = fix
        
    # CM Position Data
    X = np.load(folder + '/CMx_' + fix + '.npy')
    Y = np.load(folder + '/CMy_' + fix + '.npy')
    Z = np.load(folder + '/CMy_' + fix + '.npy')
    
    # Import Velocites
    Vx = np.load(folder + '/Vx_' + fix + '.npy')
    Vy = np.load(folder + '/Vy_' + fix + '.npy')
    Vz = np.load(folder + '/Vz_' + fix + '.npy')
    
    days = np.round(days_since_distruption(folder + '/snap_'+fix+'.h5')/t_fall,2) 

    # Calculate orbital energy
    R = np.sqrt( np.power(X,2) + np.power(Y,2) + np.power(Z,2) )
    V = np.sqrt( np.power(Vx,2) + np.power(Vy,2)+ np.power(Vz,2) )
    Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
    
    unbound_mask = np.where(Orbital > 0, 1, 0) # returns 1 on unbound gas, 0 on bound
    bound_mask = np.where(Orbital < 0, 1, 0) # returns 0 on unbound gas, 1 on bound
    
    # Import other energies
    Internal = np.load(folder + '/IE_' + fix + '.npy')
    Rad_per_m = np.load(folder + '/Rad_' + fix + '.npy')
    
    # Save
    unbound_energies = []
    unbound_energies.append( np.multiply(Orbital, unbound_mask ))
    unbound_energies.append( np.multiply(Internal, unbound_mask ))
    unbound_energies.append( np.multiply(Rad_per_m, unbound_mask ))
    unbound_energies.append( np.multiply(Mass, unbound_mask ))
        
    bound_energies = []
    bound_energies.append( np.multiply(Orbital, bound_mask ))
    bound_energies.append( np.multiply(Internal, -bound_mask )) # Minus for plot
    bound_energies.append( np.multiply(Rad_per_m, -bound_mask ))
    bound_energies.append( np.multiply(Mass, bound_mask ))
    
    return unbound_energies, bound_energies, days
    
def saver(m, fix, u, b, t):
    path = 'products/energies/' + str(m) + '/'  + str(fix)
    np.savez(path, unbound=u, bound=b, time=t)

#%%
final4 = 412
frames4 = 187
start4 = final4 - frames4

final6 = 1008
frames6 = 400
start6 = final6 - frames6

if alice:
    for i in range(frames4):
        fix4 = i + start4 + 1
        u4, b4, days4 = maker(4, fix4, alice)
        saver(4, fix4, u4, b4, days4)
        
    for i in range(frames6):    
        fix6 = i + start6 + 1
        u6, b6, days6 = maker(6, fix6, alice)
        saver(6, fix6, u6, b6, days6)
else:
    fix4 = 350
    fix6 = 820
    u4, b4, days4 = maker(4, fix4, alice)
    u6, b6, days6 = maker(6, fix6, alice)
    
#%%
def masker(x, energy_type, lower, upper):
    
    if energy_type == 'Orbital':
        i = 0
    elif energy_type == 'Internal':
        i = 1
    elif energy_type == 'Radiation':
        i = 2
    
    mask = np.logical_and(x[i]>=lower,
                                x[i]<upper)
    y = []
    for row in x:
        y.append(row[mask]) # Combine the 2 masks
    return y

# Masker
lim = 105
b6_o = masker(b6, 'Orbital', -lim, 0)
u6_o = masker(u6,  'Orbital', 0, lim)

fig, ax = plt.subplots(3,2, num=1, clear=True, tight_layout = True)
fig.suptitle('Energy type Comparison', fontsize = 15, y = 0.95)

#  Orbital ####################################################################
ax[0,0].hist( b6_o[0], bins = 30,
          weights = b6_o[-1], # mass 
          color = 'k',
          label = 'Bound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[0,0].hist(  u6_o[0], bins= 30,
          weights =  u6_o[-1], 
          color = AEK,
          label ='Unbound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[0,0].set_xlim(-110, 110)
ax[0,0].set_title('Orbital')

#  Internal ####################################################################

lim = 105
b6_i = masker(b6, 'Internal', -lim, 0)
u6_i = masker(u6,  'Internal', 0, lim)

ax[1,0].hist( b6_i[1], bins = 30,
          weights = b6_i[-1], # mass 
          color = 'k',
          label = 'Bound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[1,0].hist(  u6_i[1], bins= 30,
          weights =  u6_i[-1], 
          color = AEK,
          label ='Unbound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[1,0].set_title('Internal')
#  Radiation ##################################################################

lim = 105
b6_r = masker(b6, 'Radiation', -lim, 0)
u6_r = masker(u6,  'Radiation', 0, lim)

ax[2,0].hist( b6_r[2], bins = 30,
          weights = b6_r[-1], # mass 
          color = 'k',
          label = 'Bound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[2,0].hist(  u6_r[2], bins= 30,
          weights =  u6_r[-1], 
          color = AEK,
          label ='Unbound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[2,0].set_title('Radiation')
# 4
#  Orbital ####################################################################
lim = 25
b4_o = masker(b4, 'Orbital', -lim, 0)
u4_o = masker(u4,  'Orbital', 0, lim)

ax[0,1].hist( b4_o[0], bins = 30,
          weights = b4_o[-1], # mass 
          color = 'k',
          label = 'Bound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[0,1].hist(  u4_o[0], bins= 30,
          weights =  u4_o[-1], 
          color = AEK,
          label ='Unbound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[0,1].set_title('Orbital')

#  Internal ####################################################################

lim = 25
b4_i = masker(b4, 'Internal', -lim, 0)
u4_i = masker(u4,  'Internal', 0, lim)

ax[1,1].hist( b4_i[1], bins = 30,
          weights = b4_i[-1], # mass 
          color = 'k',
          label = 'Bound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[1,1].hist(  u4_i[1], bins= 30,
          weights =  u4_i[-1], 
          color = AEK,
          label ='Unbound',
          edgecolor = 'white',
          cumulative = False,
          log=True)

ax[1,1].set_title('Internal')
#  Radiation ##################################################################
lim = 25
b4_r = masker(b4, 'Radiation', -lim, 0)
u4_r = masker(u4,  'Radiation', 0, lim)

ax[2,1].hist( b4_r[2], bins = 30,
          weights = b4_r[-1], # mass 
          color = 'k',
          label = 'Bound',
          edgecolor = 'white',
          cumulative = False,
          log = True)

ax[2,1].hist(  u4_r[2], bins= 30,
          weights =  u4_r[-1], 
          color = AEK,
          label ='Unbound',
          edgecolor = 'white',
          cumulative = False,
          log= True)

ax[2,1].set_title('Radiation')
# Extras ######################################################################

props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
fig.text(0.15, 0.95, 't/tfb: ' + str(days6), fontsize = 10,
		    color='k', bbox = props)
fig.text(0.85, 0.95, 't/tfb: ' + str(days4), fontsize = 10,
		    color='k', bbox = props)
fig.text(0.25, 0.95, '$10^6 M_\odot$', fontsize = 15,
		    color='k')
fig.text(0.74, 0.95, '$10^4 M_\odot$', fontsize = 15,
		    color='k' )
fig.text(-0.01, 0.45, 'Mass Weighted Counts', fontsize = 15,
		    color='k', verticalalignment = 'center', rotation = 90)
fig.text(0.5, -0.01, 'Specific Energy', fontsize = 15,
		    color='k' , horizontalalignment = 'center')

#%% 
plt.figure()
fig, ax = plt.subplots(1,2, num=1, clear=True, tight_layout = True)
fig.suptitle('Energy type Comparison', fontsize = 15, y=0.92)

orb_b6, bins_b6 = np.histogram(b6_o[0], 20, weights = b6_o[-1]) 
orb_u6, bins_u6 = np.histogram(u6_o[0], 20, weights = u6_o[-1])

bins_b6 = [ (bins_b6[i] + bins_b6[i-1])/2 for i in range(1,len(bins_b6)) ] 
bins_u6 = [ (bins_u6[i] + bins_u6[i-1])/2 for i in range(1,len(bins_u6)) ] 

ie_b6, _ = np.histogram(b6_i[1], 20, weights = b6_i[-1]) 
ie_u6, _ = np.histogram(u6_i[1], 20, weights = u6_i[-1]) 

rad_b6, _ = np.histogram(b6_r[2], 20, weights = b6_r[-1]) 
rad_u6, _ = np.histogram(u6_r[2], 20, weights = u6_r[-1])

ax[0].bar(bins_b6, orb_b6, 
          width = 5, color = AEK, edgecolor = 'white', label = 'Orbital')
ax[0].bar(bins_b6, ie_b6, bottom = orb_b6,
          width = 5, color = 'k', edgecolor = 'white', label = 'Internal') 
ax[0].bar(bins_b6, rad_b6, bottom = ie_b6,
          width = 5, color = 'maroon', edgecolor = 'white', label = 'Radiation')
ax[0].legend()
#
ax[0].bar(bins_u6, orb_u6, 
          width = 5, color = AEK, edgecolor = 'white')
ax[0].bar(bins_u6, ie_u6, bottom = orb_u6,
          width = 5, color = 'k', edgecolor = 'white') 
ax[0].bar(bins_u6, rad_u6, bottom = ie_u6,
          width = 5, color = 'maroon', edgecolor = 'white') 

ax[0].set_xlim(-110, 110)

orb_b4, bins_b4 = np.histogram(b4_o[0], 20, weights = b4_o[-1]) 
orb_u4, bins_u4 = np.histogram(u4_o[0], 20, weights = u4_o[-1])

bins_b4 = [ (bins_b4[i] + bins_b4[i-1])/2 for i in range(1,len(bins_b4)) ] 
bins_u4 = [ (bins_u4[i] + bins_u4[i-1])/2 for i in range(1,len(bins_u4)) ] 

ie_b4, _ = np.histogram(b4_i[1], 20, weights = b4_i[-1]) 
ie_u4, _ = np.histogram(u4_i[1], 20, weights = u4_i[-1]) 

rad_b4, _ = np.histogram(b4_r[2], 20, weights = b4_r[-1]) 
rad_u4, _ = np.histogram(u4_r[2], 20, weights = u4_r[-1])

ax[1].bar(bins_b4, orb_b4, 
          width = 1, color = AEK, edgecolor = 'white', label = 'Orbital')
ax[1].bar(bins_b4, ie_b4, bottom = orb_b4,
          width = 1, color = 'k', edgecolor = 'white', label = 'Internal') 
ax[1].bar(bins_b4, rad_b4, bottom = ie_b4,
          width = 1, color = 'maroon', edgecolor = 'white', label = 'Radiation')
ax[1].legend()
#
ax[1].bar(bins_u4, orb_u4, 
          width = 1, color = AEK, edgecolor = 'white')
ax[1].bar(bins_u4, ie_u4, bottom = orb_u4,
          width = 1, color = 'k', edgecolor = 'white') 
ax[1].bar(bins_u4, rad_u4, bottom = ie_u4,
          width = 1, color = 'maroon', edgecolor = 'white') 

ax[1].set_xlim(-21, 21)
# Extras ######################################################################

props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
fig.text(0.15, 0.9, 't/tfb: ' + str(days6), fontsize = 10,
		    color='k', bbox = props)
fig.text(0.8, 0.9, 't/tfb: ' + str(days4), fontsize = 10,
		    color='k', bbox = props)
ax[0].set_title('$10^6 M_\odot$')
ax[1].set_title('$10^4 M_\odot$')
fig.text(0.07, 0.45, 'Mass Weighted Counts', fontsize = 15,
		    color='k', verticalalignment = 'center', rotation = 90)
fig.text(0.5, 0.07, 'Specific Energy', fontsize = 15,
		    color='k' , horizontalalignment = 'center')