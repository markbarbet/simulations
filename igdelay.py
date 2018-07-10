# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:04:47 2018

@author: Mark Barbet
"""

import cantera as ct
import numpy as np
import os
import efficiency_manipulate as em
import ig_delay as ig
import matplotlib.pyplot as plt


#Range of temperatures to run simulations for
T=np.arange(1000,1250,10)

#Pressure input *******IN ATMOSPHERES ONLY*******
P=19.3
X={'CH3OH':5.0, 'O2':15.0, 'N2':80}
mechanism=[os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem.cti',os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti']
mechanism=[os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem.cti']
efficiency_manipulate=True
reactorType='cv'

if efficiency_manipulate:
    for i in np.arange(len(mechanism)):
        gas=ct.Solution(mechanism[i])
        gas.name='igDelayRateSwap_'+mechanism[i].split('\\')[-1].rstrip('.cti')
        gas=em.efficiency_rate_swap(gas)
        import soln2cti as ctiw
        new_file=ctiw.write(gas)
        #gas2=ct.Solution(os.getcwd()+'\\pym_gas.cti')
        mechanism.append(new_file)


delays=[]
for i in np.arange(len(mechanism)):
    #gas=ct.Solution(mechanism[i])
    delays.append([])
    for j in np.arange(len(T)):
        delays[i].append(ig.ignition_delay(mechanism[i],T[j],P,X,options=reactorType))
plt.figure()
colors=['b','r','k','g','c','m']
for i in np.arange(len(mechanism)):
    plt.plot(1000.0/T,np.log10(delays[i]),color=colors[i])
    

    


