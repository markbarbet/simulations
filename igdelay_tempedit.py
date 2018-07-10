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
T=np.arange(1000,1250,50)#methanol
#T=np.arange(625,1450,100)#dme
val='all'
print(val)
#Pressure input *******IN ATMOSPHERES ONLY*******
P=19.3#methanol
#P=11.0#dme
X={'CH3OH':5.0, 'O2':15.0, 'N2':80}#methanol
#X={'ch3och3':0.3,'o2':3.0,'n2':3.76*3.0} #dme
#P=2.36#propane
#T=np.arange(1300,1500,40)#propane
#X={'C3H8':1.6,'O2':8.0,'AR':90.4}#propane
mechanism=[os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem.cti',os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti']
mechanism=[os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem.cti']
#mechanism=[os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem.cti',os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem_methanolEff.cti']

#mechanism=[os.getcwd()+'\\Mechanisms\\CH4_DME\\chem.cti']
#mechanism=[os.getcwd()+'\\Mechanisms\\Aramco2.0\\chem.cti']

efficiency_manipulate=True
reactorType='cv'

if efficiency_manipulate:
    gases=[]
    for i in np.arange(len(mechanism)):
        gas=ct.Solution(mechanism[i])
        gases.append(gas)
        gas.name='igDelayRateSwap_'+mechanism[i].split('\\')[-1].rstrip('.cti')
        gas2=em.efficiency_rate_swap(gas,[val])
        gases.append(gas2)
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
    plt.semilogy(1000.0/T,delays[i],color=colors[i])
    

    


