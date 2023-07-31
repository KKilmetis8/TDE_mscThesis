#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:19:52 2023

@author: konstantinos
"""

import numpy as np
import matplotlib.pyplot as plt
AEK = '#F1C410' # Important color
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['figure.figsize'] = [8 , 11]
plt.rcParams['axes.facecolor']='whitesmoke'

def time_averager(evolution):
    
    time_averaged = np.zeros( np.shape(evolution[0]) )
    for i in range(len(evolution)): # Loop over time
        time_averaged = np.add(time_averaged, evolution[i])
        
    time_averaged = np.multiply(time_averaged, 1/len(evolution))
    return time_averaged

outflow = []
ms = [6, 4]
for m in ms:
    folder = 'products/dai/' + str(m) + '-'
    pre = 'bins'
    quantity = '-rho' # rho or vel
    # stream = False
    
    # Load
    times1 = np.load(folder + pre + quantity + '1' + '.npy')
    times2 = np.load(folder + pre + quantity + '2' + '.npy')
    times3 = np.load(folder + pre + quantity + '3' + '.npy')
    data = [times1, times2, times3]
    
    # Tavg
    data_tavg = [ time_averager(times) for times in data]
    outflow.append(data_tavg) # Hold
#%% Plotting
# Mbh = 10**m
# Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
# apocenter = 2 * Rt * Mbh**(1/3)

# r_start = 0.2 * 2 * Rt
# r_stop = apocenter
# r_num = 100 # np.abs(x_start - x_stop)
# radii = np.linspace(r_start, r_stop, num = r_num) / apocenter
# theta_num = 7
# thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
# colarr = ['red', 'darkorange', 'gold',  'lightgreen', 'green', 'blue', 'violet']
# line = 'solid'

# if stream:
#     fig, ax = plt.subplots(1,2, tight_layout=True)
# else: 
#     fig, axs = plt.subplots(1,1, tight_layout=True)
#     ax = [axs]

# for i in range(theta_num):
#     if i == 3:
#         line = 'dashed'
#         # continue
#     label = thetas[i] * 180 / np.pi
#     label = str(np.round(label,0))
#     ax[0].plot(radii, tavg[:,i], color = colarr[i], linestyle = line,
#                label = label)
#     if stream:
#         ax[1].plot(radii, tavg_stream[:,i], color=colarr[i], linestyle = line,)
# ax[0].grid()
# ax[0].set_title('$10^' + str(m) +' M_\odot$ ' + 'Stream Included')

# if stream:
#     ax[1].grid()
#     ax[1].set_title('Stream Excluded')

# ax[0].set_xlabel(r'Radial Coordinate [$r/R_a$]')
# if quantity =='-vel':
#     ax[0].set_ylabel('Mass weighted Velocity [km/s] ')
# if quantity =='-rho':
#     ax[0].set_ylabel('Mass weighted Density [g/cm3]')
    
# ax[0].set_yscale('log')
# ax[0].tick_params(axis = 'both', which = 'both', direction='in')

# if stream:
#     ax[1].set_yscale('log')
#     ax[1].set_xlabel(r'Radial Coordinate [$R_{\odot}$]')
#     ax[1].tick_params(axis = 'both', which = 'both', direction='in')

# if quantity =='-vel':
#     fig.suptitle(r'Mass weighted, binned look at \textbf{outflow} velocity', fontsize = 16)
# if quantity =='-rho':
#     fig.suptitle(r'Mass weighted, binned look at \textbf{outflow} density', fontsize = 16)
# fig.legend(loc=7, title = r'$\theta$ bins')
# fig.tight_layout()
# fig.subplots_adjust(right=0.88)
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
# fig.text(0.8, 1.05, 'This uses only 2 snapshots! \n It will be deployed at scale.',
#           transform=ax[0].transAxes, bbox = props)

#%%
Mbh = 10**m
Rt =  Mbh**(1/3) # Msol = 1, Rsol = 1
apocenter = 2 * Rt * Mbh**(1/3)

r_start = 0.2 * 2 * Rt
r_stop = 4 * apocenter
r_num = 400 # np.abs(x_start - x_stop)
radii = np.linspace(r_start, r_stop, num = r_num) / apocenter
theta_num = 7
thetas = np.linspace(-np.pi/2, np.pi/2, num = theta_num)
colarr = ['red', 'darkorange', AEK,  'k', 'green', 'blue', 'violet']
line = 'solid'
m6 = r'$10^6 M_\odot$'
m4 = r'$10^4 M_\odot$'

titles4 = ['0.5 - 0.71', '0.71 - 1', '1 - 1.3'] 
titles6 = ['0.5 - 0.91', '0.91 - 1.39', '1.39 - 1.61']
labels = []
fig, axs = plt.subplots(3,2, tight_layout=True, sharex=True, sharey=True)

# Ax lims and scale
if quantity == '-rho':
    #axs[0,0].set_ylim(1e0,6e4)
    axs[0,0].set_yscale('log')
if quantity == '-vel':
    axs[0,0].set_ylim(1e0,6e4)
    axs[0,0].set_yscale('log')

for k in range(3):
    plt.figure()
    line = 'solid'
    for i in range(theta_num):
        if i == 4:
            line = 'dashed'
            # continue
        
        # Make labels
        if k == 0:
            label = thetas[i] * 180 / np.pi
            labels.append(str(np.round(label,0)))
            
            axs[k,0].plot(radii, outflow[0][k][:,i], color = colarr[i], linestyle = line, 
                          label = labels[i])
            axs[k,1].plot(radii, outflow[1][k][:,i], color = colarr[i], linestyle = line)
        else:    
            axs[k,0].plot(radii, outflow[0][k][:,i], color = colarr[i], linestyle = line)
            axs[k,1].plot(radii, outflow[1][k][:,i], color = colarr[i], linestyle = line) 
            
        # Make pretty
        axs[k,0].tick_params(axis = 'both', which = 'both', direction='in')
        axs[k,1].tick_params(axis = 'both', which = 'both', direction='in')
        axs[k,0].grid()
        axs[k,1].grid()
        axs[k,0].set_title(m6 + ' for ' + titles6[k] + r' $t_{FB}$', fontsize = 18)
        axs[k,1].set_title(m4 + ' for ' + titles4[k] + r' $t_{FB}$', fontsize = 18)
        axs[k,0].yaxis.set_tick_params(labelsize=13)
        axs[k,1].yaxis.set_tick_params(labelsize=13)
        axs[2,0].xaxis.set_tick_params(labelsize=15)
        axs[2,1].xaxis.set_tick_params(labelsize=15)
  


# Big M
# txt_x = 0.01
# txt_y = 1.03 
# fig.text(txt_x, txt_y, '$10^6 M_\odot$', 
#         fontsize = 30,
#         fontfamily = 'bold',
#         color = 'k',
#         transform = axs[0,0].transAxes)
# fig.text(txt_x + 1.78 , txt_y, '$10^4 M_\odot$', 
#         fontsize = 30,
#         fontfamily = 'bold',
#         color = 'k',
#         transform= axs[0,0].transAxes)


fig.legend(loc = 'lower center', ncols = 7, columnspacing=0.8, fontsize = 14,
           bbox_to_anchor = (0.02,-0.038,1,1))
fig.supxlabel('Radial Coordinate [$r/R_a$]', fontsize = 20)
if quantity =='-vel':
    fig.suptitle('Outflow Velocity', 
              fontsize =  25)
    fig.supylabel('Mass weighted Velocity [km/s]', 
              fontsize =  20)
if quantity =='-rho':
    fig.suptitle('Outflow Density', 
              fontsize =  25)
    fig.supylabel('Mass weighted Density [g/cm$^3$] ', 
              fontsize =  25)