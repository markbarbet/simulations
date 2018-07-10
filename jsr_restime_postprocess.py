# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 16:49:22 2017

@author: Mark Barbet

Code to do post-processing of reactor restime vs Temperature data
"""

import os
import matplotlib.pyplot as plt
tempresults=results
import numpy as np
mech=[]
for i in tempresults:
    if i.mechanism not in mech:
        mech.append(i.mechanism)

for j in np.arange(len(mech)):
    tempRes=[]
    tempTemp=[]
    plt.figure()
    for i in tempresults:
        if i.mechanism==mech[j]:
            tempTemp.append(i.solution['temperature'])
            tempRes.append(i.residence_time)
            currentFuel=i.fuel
            currentPhi=str(i.phi)
    plt.semilogx(tempRes,tempTemp)
    plt.xlabel('Residence Time')
    plt.ylabel('Temperature (K)')
        
    plt.title('Temperature vs. Residence Time for '+mech[j].split('\\')[-1]+'_ethene')
    #plt.legend(bbox_to_anchor=(0.5,-0.2),loc='upper center',frameon=True,prop={'size':10})
    #from matplotlib.font_manager import FontProperties
    #fontP = FontProperties()
    #fontP.set_size('small')
    #plt.tight_layout(pad=1.5)
    #plt.legend([pl], 'Flame Speeds for '+j, prop = fontP)
    
    #plt.savefig(os.getcwd()+'\\figures\\'+currentFuel+currentPhi+'_'+mech[j].split('\\')[-1]+'_resTemps.png',dpi=1200,bbox_inches='tight')
    if 'appendedinput' not in mech[j]:
        baseTemp=tempTemp
        baseRes=tempRes
    if 'appendedinput' in mech[j]:
        plt.figure()
        plt.semilogx(np.array(tempRes),np.array(tempTemp)-np.array(baseTemp))
        plt.figure()
        plt.semilogx(baseRes,baseTemp,tempRes,tempTemp)
        plt.title('Temperature vs. Residence Time for '+mech[j].split('\\')[-1]+'_'+fuel)
        plt.xlabel('Residence Time')
        plt.ylabel('Temperature (K)')
        plt.savefig(os.getcwd()+'\\figures\\'+currentFuel+currentPhi+'_'+mech[j].split('\\')[-1]+'_resTemps.png',dpi=1200,bbox_inches='tight')
