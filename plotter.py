# -*- coding: utf-8 -*-
"""
Created on Tue Jun 05 15:10:54 2018

@author: HP USER
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def plotSpecies(results,species,conditions):
    
    temps=[]
    for i in results:
        if i.solution['temperature'].values[0] not in temps:
            temps.append(i.solution['temperature'].values[0])
            
    filename=[]
    for i in results:
        if i.mechanism not in filename:
            filename.append(i.mechanism)
    res=[]
    for i in results:
        if i.residence_time not in res:
            res.append(i.residence_time)
            
            
    length=len(conditions)
    length2=len(res)
    
    
    sortedresults=[results[i::length] for i in range(length)]
    
    for j in np.arange(len(sortedresults)):
        sortedresults[j]=[sortedresults[j][i::length2] for i in range(length2)]
        
    for j in np.arange(len(sortedresults)):
        for f in np.arange(len(sortedresults[j])):
            #sortedresults[j][f]=[sortedresults[j][f][i::len(species)] for i in range(len(species))]
            sortedresults[j][f]=[sortedresults[j][f][i::len(filename)] for i in range(len(filename))]
#    for j in np.arange(len(sortedresults)):
#        for f in np.arange(len(sortedresults[j])):
#            for poop in np.arange(len(sortedresults[j][f])):
#                sortedresults[j][f][poop]=[sortedresults[j][f][poop][i::len(filename)] for i in range(len(filename))]
#                
    #return sortedresults       
    for j in np.arange(len(sortedresults)):
        for f in np.arange(len(sortedresults[j])):
            for poop in np.arange(len(species)):
                plt.figure()
                for c in np.arange(len(filename)):      
                    specieslist=[]
                    for d in np.arange(len(temps)):
                        #print(j,f,poop,c,d)
                        specieslist.append(sortedresults[j][f][c][d].solution[species[poop]].values[0])
                    plt.plot(temps,specieslist)
                    plt.xlabel('Temperature (K)')
                    plt.ylabel('Mole Fraction')
                    plt.title(str(species[poop]+' with residence time: '+str(res[f])+' sec'))
                plt.savefig('speciation'+species[poop]+'_conditions'+str(j)+'_resTime'+str(f)+'.pdf',dpi=1200,bbox_inches='tight')
            
    
    return sortedresults
#    for c in np.arange(len(sortedresults)):
#        for s in np.arange(len(species)):
#            
#            for r in np.arange(len(res)):
#                plt.figure()
#                plt.plot(temps,sortedresults[c][r])
    
def plotsens(sortedresults,observables,conditions,filename,temps,top10):
    
    linetypes=['k-','b-','r-','g-','c-','m-','y-','k-.','b-.','r-.','g-.','c-.','m-.','y-.','k--','b--','r--','g--','c--','m--','y--',
               'k:','b:','r:','g:','c:','m:','y:']
    count=0
    for j in np.arange(len(sortedresults)):
        for f in np.arange(len(sortedresults[j])):
            for o in np.arange(len(observables)):
                fig=plt.figure()
                ax=plt.subplot(111)
                ax.set_xlabel('Temperature (K)')
                ax.set_ylabel(r'Sensitivity coefficient, $\frac{\partial\mathrm{ln}S_u^0}{\partial\mathrm{ln}k}$')
                plt.title('Top '+str(top10)+' sensitivities for observable '+observables[o])
                for c in np.arange(len(filename)):
                    
                        sens=[]
                        sens=pd.DataFrame(sortedresults[j][f][c][0].Index[1],columns=['rxn'])
                        for d in np.arange(len(temps)):
                            #print(j,f,c,d,o)
                            sens=pd.concat([sens, pd.DataFrame(sortedresults[j][f][c][d].k_sens[0],columns=[temps[d]]*len(observables)).iloc[:,o]], axis=1)
                        maxes=sens.loc[:, sens.columns != 'rxn'].max(axis=1)
                        maxes=maxes.nlargest(top10)
                        max_index=np.array(maxes.index)
                        
                        for i in max_index:
                            plt.plot(temps,sens.iloc[i][1:],linetypes[count],label=sens['rxn'][i]+', '+os.path.splitext(filename[c])[0].split('\\')[-1])
                            count=count+1
                        box = ax.get_position()
                        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                                  fancybox=False, shadow=False, ncol=1)
                plt.savefig('sensitivities_conditions'+str(j)+'_resTime'+str(f)+'_'+observables[o]+'.pdf',dpi=1200,bbox_inches='tight')
