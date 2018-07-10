# -*- coding: utf-8 -*-
"""
Created on Thu Nov 02 11:46:43 2017

@author: Mark Barbet
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import copy
import os
y=results


mechlist=[]
for i in y:
    if i.mechanism not in mechlist:
        mechlist.append(i.mechanism)
flamespeed=np.zeros(np.shape(np.array(mechlist)))    
for i in np.arange(len(mechlist)):
    for j in y:
        if j.mechanism==mechlist[i] and j.solution['u'][0]>flamespeed[i]:
            flamespeed[i]=j.solution['u'][0]
            
tempmech=copy.deepcopy(mechlist)
#tempflmspeed=flamespeed
order=np.argsort(flamespeed)
for i in np.arange(len(order)):
    mechlist[i]=tempmech[order[i]]
    
markerlist=['k-','r:','g-','b-','k-.','r--','g-.','k:','r-','g--','k--','r-.','g:',]

#colors=['k','g','r','b','m','c','y']
#styles=[':','--','-','-.']
#markerlist=[]
#for i in styles:
#    for j in colors:
#        markerlist.append(j+i)
phiList=[]
for i in y:
    if i.phi not in phiList:
        phiList.append(i.phi)
reducedMechs=[]
for i in mechlist:
    reducedMechs.append(i.split('\\')[-1].rstrip('.cti'))
flamespeedResults=pd.DataFrame(columns=reducedMechs+['phi'])
flamespeedResults['phi']=phiList
for i in mechlist:
    tempflmspeed=[]
    #phis=[]
    for j in y:
        for p in phiList:
            if j.mechanism==i and j.phi==p:
                tempflmspeed.append(j.solution['u'][0])
                #phis.append(j.phi)
    flamespeedResults[i.split('\\')[-1].rstrip('.cti')]=tempflmspeed
plt.figure()
    
for i in reducedMechs:
    #count=count+1
    
    
    plt.plot(flamespeedResults['phi'],flamespeedResults[i],markerlist[reducedMechs.index(i)],linewidth=2,label=str(reducedMechs.index(i)+1)+'. '+i.split('\\')[-1].rstrip('.cti'))
    plt.xlabel('Phi')
    plt.ylabel('Flame Speed (m/s)')
        
    plt.title('Flame Speeds for '+y[0].fuel)
    plt.legend(bbox_to_anchor=(0.5,-0.2),loc='upper center',frameon=True,prop={'size':10})
    #from matplotlib.font_manager import FontProperties
    #fontP = FontProperties()
    #fontP.set_size('small')
    plt.tight_layout(pad=1.5)
    #plt.legend([pl], 'Flame Speeds for '+j, prop = fontP)
    
    plt.savefig(os.getcwd()+'\\figures\\'+y[0].fuel+'\\'+y[0].fuel+'_flamespeed.png',dpi=1200,bbox_inches='tight')
    
plt.figure()
diffmechs=[]
for i in mechlist:
    if 'appendedinput' in i:
        diffmechs.append(i)
    if 'appendedinput' not in i:
        basemech=i

for i in diffmechs:
    plt.plot(flamespeedResults['phi'],np.divide((flamespeedResults[i.split('\\')[-1].rstrip('.cti')]-flamespeedResults[basemech.split('\\')[-1].rstrip('.cti')]),flamespeedResults[basemech.split('\\')[-1].rstrip('.cti')]),markerlist[reducedMechs.index(i.split('\\')[-1].rstrip('.cti'))],linewidth=2,label=str(mechlist.index(i)+1)+'. '+i.split('\\')[-1].rstrip('.cti'))
    plt.xlabel('Phi')
    plt.ylabel('Flame Speed Relative Deviation')
        
    plt.title('Flame Speed Relative Difference for '+y[0].fuel)
    plt.legend(bbox_to_anchor=(0.5,-0.2),loc='upper center',frameon=True,prop={'size':10})
    #from matplotlib.font_manager import FontProperties
    #fontP = FontProperties()
    #fontP.set_size('small')
    #plt.tight_layout(pad=1.5)
    #plt.legend([pl], 'Flame Speeds for '+j, prop = fontP)
    
    plt.savefig(os.getcwd()+'\\figures\\'+y[0].fuel+'\\'+y[0].fuel+'_flamespeed_rel_diffs.png',dpi=300,bbox_inches='tight')

