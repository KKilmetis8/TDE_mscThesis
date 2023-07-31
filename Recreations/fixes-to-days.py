#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:41:32 2023

@author: konstantinos
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] = [6 , 6]
plt.rcParams['font.family'] = 'Serif'
from extractors.time_extractor import days_since_distruption
##
days = []
fixes = [820,881,925]
for fix in fixes:
    # Progress Check
    days.append(str(np.round(days_since_distruption(fix+'/snap_'+fix+'.h5'),1)))

plt.plot(fixes, days, c='k')
plt.xlabel('Fixes', fontsize = 14)
plt.ylabel('Days', fontsize = 14)
plt.title('Conversion Chart: Fixes to Days', fontsize = 18)
plt.grid()
plt.savefig('fix-to-days.png')
        