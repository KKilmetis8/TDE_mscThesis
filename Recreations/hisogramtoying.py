# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 15:19:05 2022

@author: Konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
from src.Extractors.time_extractor import linear_fit_days

plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 300
# Plot size
plt.rcParams['axes.facecolor']='whitesmoke'
plt.rcParams['figure.figsize'] = [3.0, 3.0]
plt.rcParams['font.family'] = 'Times New Roman'
#%%
# Conversion constants
G = 6.6743e-11 # SI
Msol = 1.98847e30 # kg
Rsol = 6.957e8 # m
t = np.sqrt(Rsol**3 / (Msol*G )) # Follows from G=1
Mbh = 1e6 # * Msol
# Need these for the PW potential
c = 3e8 * t/Rsol # c in simulator units.
rg = 2*Mbh/c**2
t_fall = 40 * (Mbh/1e6)**(0.5) # days EMR+20 p13

# Data Load
fix = '820'
fixdays = str(linear_fit_days(int(fix)))
mass = np.load(fix + '/Mass_' + fix + '.npy')
X = np.load(fix + '/CMx_' + fix + '.npy')
Y = np.load(fix + '/CMy_' + fix + '.npy')
Z = np.load(fix + '/CMz_' + fix + '.npy')
Vx = np.load(fix + '/Vx_' + fix + '.npy')
Vy = np.load(fix + '/Vy_' + fix + '.npy')
Vz = np.load(fix + '/Vz_' + fix + '.npy')
R = np.sqrt( np.power(X,2) + np.power(Y,2)+ np.power(Z,2))
V = np.sqrt( np.power(Vx,2) + np.power(Vy,2)+ np.power(Vz,2))
Orbital = (0.5 * V**2 ) - Mbh / (R-rg)
Bound = np.where(Orbital < 0, 1, 0)
# Spec. energy for bound
spec_bound =  Orbital*Bound # *mass # The minus is me cheating

# Get unbound -1 turns 0 -> -1 and 1->0. 
# For a weird potential reason, 1008 gets negative energies
Unbound = (Bound - 1) * (-1)

# Spec. energy for unbound
spec_unbound = Orbital * Unbound # * mass

# Mask extreme values.
masking = True
if masking:
    lower_bound = -150 # 1000 for specific, 0.75e-5 for simunits
    upper_bound = 150
    bound_mask = np.logical_and(spec_bound>=lower_bound,
                                spec_bound<upper_bound)
    unbound_mask = np.logical_and(spec_unbound>=lower_bound,
                                  spec_unbound<upper_bound)
    mask = np.logical_and(bound_mask, unbound_mask) # Combine the 2 masks
    
    mass = mass[mask]
    Bound = Bound[mask]
    Unbound = Unbound[mask]
    spec_bound = spec_bound[mask]
    spec_unbound = spec_unbound[mask]
#%% Histograms

# Calculate bound and unbound mass
mass_bound = Bound*mass
mass_unbound = np.abs(Unbound*mass)

# Plot
# plt.hist( (spec_bound, spec_unbound), bins=50,
#           weights = (mass_bound, mass_unbound), 
#           color = ('mediumslateblue', 'tomato'),
#           label = ('Bound', 'Unbound'),
#           cumulative=True,
#           log=True)
dm, de, _ = plt.hist( spec_bound, bins=50, # BEEG PAPATZILIKI
          weights = mass_bound, 
          color = '#001158',
          label = 'Bound',
          cumulative = False,
          log=True,
          alpha = 0.9)
plt.hist(  spec_unbound, bins=50,
          weights =  mass_unbound, 
          color = '#be1908',
          label ='Unbound',
          cumulative = False,
          log=True,
          alpha = 0.9)

# Plotting boilerplate
plt.legend(loc = 'lower center')
# plt.xlim(-250, 250) 
# plt.ylim(1e-3,1e-1)
plt.title('Mass Weighted Histogram of Specific energies after ' 
          + fixdays + ' days', pad=10, fontsize=10)
plt.xlabel(r'Specific Energy', fontsize=10)
plt.ylabel(r'Mass $\cdot$ Counts $\left[M_{\odot} \cdot \# \right]$', 
           fontsize=10)

# Derivative 
np.save('dm'+fix, dm)

#%% Energy Bars
 
# -1 turns 0 -> -1 and 1->0. *-1 gets us positive
Unbound =  (Bound - 1) * (-1)


# Potential energies are negative, therefore there is a minus there
pot_b = -np.dot( pot*mass, Bound)  
pot_u = -np.dot( pot*mass, Unbound) 
print( '%.1e' % pot_u)
kin_b = np.dot(kin*mass, Bound) 
kin_u = np.dot(kin*mass, Unbound) 

plt.bar( x = ('Bound Pot', 'Unbound Pot', 'Bound Kin', 'Unbound Kin'),
         height = (pot_b,pot_u, kin_b, kin_u),
         width = 1,
         color =('mediumslateblue', 'tomato', 'darkslateblue', 'maroon'),
         log=False
         )
plt.xticks(fontsize=8, rotation=15)
#plt.ylim(-1e43, 1e43)
plt.ylabel('Specific Energy [sim units]')
plt.title('Potential \& Kinetic specific energies after ' 
          + fixdays + ' days', pad=10, fontsize=15)

# Total Energy
plt.figure()
plt.title('Total Energy in bound and unbound parts after ' 
          + fixdays + ' days', pad=10, fontsize=15)
plt.bar( x = ('Bound', 'Unbound'),
         height = (pot_b , pot_u ),
         color = ('mediumslateblue', 'tomato'),
         label = 'Potential'
         )
plt.bar( x = ('Bound', 'Unbound'),
         height = ( kin_b, kin_u),
         color = ('darkslateblue', 'maroon'),
         bottom = (pot_b, pot_u),
         label = 'Kinetic'
         )
plt.legend()

#%% Figure 2

# Data in.
fixes = ['820', '881' , '925', '950']#, '1008']
fixdays = [ str(linear_fit_days(int(fix))) for fix in fixes ]
bounds = [ np.load(fix+'/BoundCM_'+fix+'.npy') for fix in fixes]
masses = [ np.load(fix+'/Mass_'+fix+'.npy') for fix in fixes]


#%%
def get_bound_mass(Mass, Bound):
    bound_mass = np.dot(Mass,Bound)
    return bound_mass

bound_masses = [get_bound_mass(masses[i], bounds[i])
                for i in range(len(fixes)) ]

# Plot size
plt.rcParams['figure.figsize'] = [4.0, 4.0]
plt.figure()
plt.plot( fixdays , bound_masses, 'o-',color = 'purple')
plt.title('Bound Mass')
plt.xlabel('Days since distruption')
plt.ylabel(r'Bound Mass $\left[ M_{\odot} \right]$')
