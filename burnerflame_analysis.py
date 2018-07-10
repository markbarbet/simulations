# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 13:26:06 2017

@author: Mark Barbet
"""

#File for BurnerFlame Analysis
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

observables=obs
tempresults=results
mechList=[]
#sens=

for i in tempresults:
    if i.mechanism not in mechList:
        mechList.append(i.mechanism)
markerlist=['k-','k:','k-.','k--','r-','r:','r-.','r--']        
for j in observables:
    print(j)
    for i in np.arange(len(tempresults)):
        if i==0:
            baseZ=np.array(tempresults[i].Index[0])
            baseObs=np.array(tempresults[i].solution[j])
        else:
            plt.figure()
            plt.plot(baseZ,baseObs,tempresults[i].Index[0],tempresults[i].solution[j])
            #plt.title('Mole Fractions for '+i.fuel+' at phi='+str(i.phi))
            plt.title('Mole Fractions for '+j+' in '+tempresults[i].mechanism.split('\\')[-2].split('.')[0]+'_'+tempresults[i].mechanism.split('\\')[-1].split('.')[0])
            plt.xlabel('z (m)')
            plt.ylabel(j)
            #plt.legend(loc='bottom right',frameon=True,prop={'size':10})
            plt.savefig(os.getcwd()+'\\figures\\'+tempresults[i].mechanism.split('\\')[-2].split('.')[0]+'_'+tempresults[i].mechanism.split('\\')[-1].split('.')[0]+'_'+tempresults[i].fuel+str(tempresults[i].phi)+'_'+j+'.png'
            ,dpi=1200,bbox_inches="tight",format='png')
for j in np.arange(len(observables)):
    plt.figure()
    plt.title('Mole Fractions for ' +observables[j])
    plt.xlabel('z (m)')
    plt.ylabel(observables[j])
    
    for i in np.arange(len(tempresults)):
        
        plt.plot(tempresults[i].Index[0],tempresults[i].solution[observables[j]],markerlist[i],label=tempresults[i].mechanism.split('\\')[-1])
        
    plt.savefig(os.getcwd()+'\\figures\\'+'burnerFlame_konnov_'+observables[j]+'.png',dpi=1200,bbox_inches='tight',format='png')
    plt.legend(loc=4,frameon=True,prop={'size':10})    
if sens==1:
    tempsens=[]
    sensitivities=[]
    for j in tempresults:
        sensitivities.append(pd.DataFrame(columns=tempresults[0].Index[2],index=j.Index[1]))
    for i in np.arange(len(tempresults)):
        
        tempsens=tempresults[i].x_slice_ksens(tempresults[i].Index[0][-1])
        sensitivities[i]=pd.DataFrame(tempsens,index=tempresults[i].Index[1])
        sensitivities[i].columns=tempresults[i].Index[2]
        
        
        sortbool=[False]*len(tempresults[0].Index[1])
        original_length=len(tempresults[0].Index[1])
    #for i in tempresults:
        
    for i in np.arange(len(tempresults)):
        if len(tempresults[i].Index[1])> original_length:
            sortbool=[False]*len(tempresults[i].Index[1])
        for k in tempresults[i].Index[2]:
            if len(sortbool)>len(tempresults[0].Index[1]):
                for j in np.arange(len(sensitivities[i])):
                    if tempresults[i].Index[1][j] in tempresults[0].Index[1]:
                        sortbool[j]=False
                    else:
                        sortbool[j]=True
                termolecularSubset=pd.DataFrame(sensitivities[i][k][sortbool])
                termolecularSubset['abs']=np.abs(termolecularSubset[k])
                temp=pd.DataFrame(termolecularSubset.nlargest(10,'abs'))
                del(temp['abs'])
                try:
                    temp.plot.barh(title="Sensitivities for "+k+'_'+tempresults[i].mechanism.split('\\')[-1].split('.')[0]+'_'+tempresults[i].fuel,legend=None)
                    #threshold = 0.001
                    #firstColumn = termolecularSubset.columns[0]
                    
                    # For plotting, collect only those steps that are above the threshold
                    # Otherwise, the y-axis gets crowded and illegible
                    #sensitivitiesSubset =termolecularSubset[termolecularSubset[firstColumn].abs() > threshold]
                    #print(sensitivitiesSubset)
                    #indicesMeetingThreshold =sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False).index
                    #sensitivitiesSubset.loc[indicesMeetingThreshold].plot.barh(title="Sensitivities for "+i.mechanism.split('\\')[-1].split('.')[0]+'_'+i.fuel+" at phi="+str(i.phi),
                    #                                                      legend=None)
                    plt.gca().invert_yaxis()
                
                    plt.rcParams.update({'axes.labelsize': 20})
                    plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');
                
                    # Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
                    plt.savefig(os.getcwd()+'\\figures\\'+tempresults[i].mechanism.split('\\')[-1].split('.')[0]+'_'+tempresults[i].fuel+'_'+k+'.png'
                    ,dpi=1200,bbox_inches="tight",format='png')
                except:
                    print('Problem plotting for '+k+' with mechanism '+tempresults[i].mechanism.split('\\')[-1].split('.')[0])
        
        
        
        