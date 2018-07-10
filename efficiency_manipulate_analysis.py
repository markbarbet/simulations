# -*- coding: utf-8 -*-
"""
Created on Mon Jan 08 17:25:20 2018

@author: Mark Barbet
"""

#Analysis of results for flamespeeds from efficiency script

import numpy as np
import pandas as pd
import os
import copy
import matplotlib.pyplot as plt


datacopy=copy.deepcopy(results)
baseline_flmspd=[]
modified_flmspd=[]
phi=[]
phi2=[]
for i in np.arange(len(datacopy)):
    try:
        if i%2==0:
            modified_flmspd.append(datacopy[i].solution['u'][0])
            phi.append(datacopy[i].phi)
        elif i%2!=0:
            baseline_flmspd.append(datacopy[i].solution['u'][0])
            phi2.append(datacopy[i].phi)
    except:
        
        print('Index ' +str(i)+' failed.')
plt.figure()
plt.plot(phi2,baseline_flmspd,phi,modified_flmspd)
    