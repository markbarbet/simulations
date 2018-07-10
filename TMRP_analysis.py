# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 11:51:01 2017

@author: Mark Barbet
"""

"""
This file is for results processing from the free flame model data object.  
Results are created from the appending of model data objects to a list of results
Use this for plotting and analysis.  This file assumes the existance of a results object in the workspace.
"""
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
graphs=False
a=[]
resultsSubset=[]
for i in results:
    if i!='failed':
        resultsSubset.append(i)
for d in resultsSubset:
    a.append(d)
phi_list=[]
for i in np.arange(len(resultsSubset)):
    phi_list.append(resultsSubset[i].phi)
b=pd.DataFrame(data=phi_list)
b.columns=['Phi']
u=[]
for i in np.arange(len(resultsSubset)):
    u.append(resultsSubset[i].solution['u'][0])
b['u']=u
#b['solution']=results.values()

mech=[]

for j in np.arange(len(resultsSubset)):
    mech.append(resultsSubset[j].mechanism)
mechList=[]
for i in mech:
    if i not in mechList:
        mechList.append(i)
solutions=[]
for i in np.arange(len(resultsSubset)):
    solutions.append(resultsSubset[i])
b['mechanism']=mech
b['solution']=solutions
b=b.sort_values(by='Phi', ascending=1)
for k in mechList:
    subset=b[b['mechanism']==k]
    #plt.plot(subset['Phi'],subset['u'])
for j in resultsSubset:
    if 'dict' in str(type(j.fuel)):
        tempstr=''
        for i in j.fuel.keys():
            tempstr=tempstr+i+'+'
        tempstr=tempstr.rstrip('+')
        j.fuel=tempstr
for j in np.arange(len(fuels)):
    if 'dict' in str(type(fuels[j])):
        tempstr=''
        for i in fuels[j].keys():
            tempstr=tempstr+i+'+'
        tempstr=tempstr.rstrip('+')
        fuels[j]=tempstr
for i in fuels:
    b[i]=[False]*len(resultsSubset)
for i in np.arange(len(resultsSubset)):
    for k in fuels:
        if b['solution'][i].fuel==k:
            b[k][i]=True
markerlist=['k-','r:','g-','b-','k-.','r--','g-.','k:','r-','g--','k--','r-.','g:','b--','c-']

for j in fuels:
    plt.figure()
    for k in mechList:
        subset=b[b['mechanism']==k]
        subset=subset[subset[j]==True]
        plt.plot(subset['Phi'],subset['u'],markerlist[mechList.index(k)],linewidth=2,label=k.split('\\')[-1])
        plt.xlabel('Phi')
        plt.ylabel('Flame Speed (m/s)')
        
        plt.title('Flame Speeds for '+j)
    plt.legend(bbox_to_anchor=(0.5,-0.2),loc='upper center',frameon=True,prop={'size':10})
    #from matplotlib.font_manager import FontProperties
    #fontP = FontProperties()
    #fontP.set_size('small')
    plt.tight_layout(pad=1.5)
    #plt.legend([pl], 'Flame Speeds for '+j, prop = fontP)
    
    #plt.savefig(os.getcwd()+'\\figures\\'+fuels[0]+'\\'+j+'_uscmech_flamespeed.png',dpi=1200,bbox_inches='tight')
    plt.savefig(os.getcwd()+'\\figures\\methanol_colliders.pdf',dpi=1200,bbox_inches='tight')
    
y=[]

for g in resultsSubset:
    #if g.phi==1:
        y.append(g)

    
sortbool=[False]*len(y[0].flamespeed_sens)
lenlist=[]
for i in y:
    if len(i.flamespeed_sens) not in lenlist:
        lenlist.append(len(i.flamespeed_sens))
minLength=min(lenlist)
#for i in y:
#    if len(i.flamespeed_sens)> len(sortbool):
#        sortbool=[False]*len(i.flamespeed_sens)
top30=[]    
for i in y:
    sortbool=[False]*len(i.flamespeed_sens)
    if len(sortbool)>minLength:
    #if len(sortbool)==len(i.flamespeed_sens):
        for j in np.arange(len(i.flamespeed_sens)):
            if i.flamespeed_sens.index[j] in y[0].flamespeed_sens.index:
                sortbool[j]=False
            else:
                sortbool[j]=True
        termolecularSubset=pd.DataFrame(i.flamespeed_sens[sortbool])
        termolecularSubset['abs']=np.abs(termolecularSubset['Su'])
        temp=pd.DataFrame(termolecularSubset.nlargest(20,'abs'))
        del(temp['abs'])
        temp.plot.barh(title="Sensitivities for "+i.mechanism.split('\\')[-1].split('.')[0]+'_'+i.fuel+" at phi="+str(i.phi),legend=None)
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
        #plt.savefig(os.getcwd()+'\\figures\\'+i.fuel+'\\'+i.mechanism.split('\\')[-1].split('.')[0]+'_'+i.fuel+'_phi'+str(i.phi).replace('.','')+'.png'
       # ,dpi=300,bbox_inches="tight",format='png')
        for p in np.arange(len(temp)):
            if temp.index[p] not in top30:
                top30.append(temp.index[p])

phiList=[]
for i in y:
    if i.phi not in phiList:
        phiList.append(i.phi)
radicals=['OH','H','O','HO2']
appendedmechs=[]
for i in mechList:
    if 'appendedinput' in i:
        appendedmechs.append(i)
forwardlist=[]
reverselist=[]
for i in appendedmechs:
    if 'forward' in i:
        forwardlist.append(i)
    elif 'forward' not in i:
        reverselist.append(i)
if 'forward_rates' in dir(y[0]) and graphs==True:
    
    for i in mechList:
        rates=[]
        for j in y:
            if j.mechanism==i and j.phi in [min(phiList),max(phiList),1.0]:
                #for x in np.arange(10):
                #for x in np.arange(len(reference.Index[1])):
                #plt.figure()
                for reaction in np.arange(len(j.Index[1])):
                    if j.Index[1][reaction] in top30:
                        fig = plt.figure()
                        ax1 = fig.add_subplot(111)
                        ax1.plot(j.Index[0],j.net_rates_of_progress[reaction,:])
                        ax1.set_ylabel('w')
                        
                        ax2 = ax1.twinx()
                        ax2.plot(j.Index[0],j.solution['T'], 'r-')
                        ax2.set_ylabel('T', color='r')
                        plt.title(j.Index[1][reaction]+' for '+fuels[0]+'_phi='+str(j.phi)+' in '+i.split('\\')[-1].rstrip('.cti'))
                        plt.xlim([0.102,0.115])
                        #plt.savefig(os.getcwd()+'\\figures\\'+fuels[0]+'\\'+j.mechanism.split('\\')[-1].rstrip('.cti')+'_'+'_phi'+str(j.phi).replace('.','')+j.Index[1][reaction].replace('=','-').replace('>','-').replace('<','-')+'.png'
                        #                ,dpi=300,bbox_inches="tight",format='png')
                        plt.figure()
                        plt.plot(j.solution['T'],j.net_rates_of_progress[reaction,:])
                        plt.xlabel('T (K)')
                        plt.ylabel('w')
                        plt.title(j.Index[1][reaction]+' for '+fuels[0]+'_phi='+str(j.phi)+' in '+i.split('\\')[-1].rstrip('.cti'))
                        #plt.xlim([0.102,0.115])
                        #plt.savefig(os.getcwd()+'\\figures\\'+fuels[0]+'\\'+j.mechanism.split('\\')[-1].rstrip('.cti')+'_'+'_phi'+str(j.phi).replace('.','')+'_tempProf_'+j.Index[1][reaction].replace('=','-').replace('>','-').replace('<','-')+'.png'
                                    #    ,dpi=300,bbox_inches="tight",format='png')
#                        for tl in ax2.get_yticklabels():
#                            tl.set_color('r')
                        #plt.plot(reference.Index[0],rates[-1][x,:])
                        #plt.title(reference.Index[1][x])
#                for species in np.arange(len(radicals)):
#                    fig=plt.figure()
#                    if 'appendedinput' in j.mechanism:
#                        for t in y:
#                            if 'appendedinput' not in t.mechanism and t.phi==j.phi:
#                                baseX=t.solution[radicals[species]]
#                                baseZ=t.Index[0]
#                                baseT=t.solution['T']
#                        ax1=fig.add_subplot(111)
#                        ax1.plot(baseZ,baseX,j.Index[0],j.solution[radicals[species]])
#                        #ax1.xlabel('z')
#                        ax1.set_ylabel('X')
#                        ax2 = ax1.twinx()
#                        ax2.plot(baseZ,baseT,'r--')
#                        ax2.plot(j.Index[0],j.solution['T'],'r-.')
#                        ax2.set_ylabel('T',color='r')
#                        ax2.set_ylim([2450,2600])
#                        plt.title(radicals[species]+' for '+fuels[0]+'_phi='+str(j.phi)+' in '+i.split('\\')[-1].rstrip('.cti'))
#                        plt.savefig(os.getcwd()+'\\figures\\'+fuels[0]+'\\'+j.mechanism.split('\\')[-1].rstrip('.cti')+'_'+'_phi'+str(j.phi).replace('.','')+radicals[species]+'_speciesProfile'+'.png'
#                                        ,dpi=300,bbox_inches="tight",format='png')
#for species in np.arange(len(radicals)):
#    for m in forwardlist:
#        for x in [min(phiList),max(phiList),1.0]: 
#            fig=plt.figure()
#            for t in y:
#                if t.phi==x and m==t.mechanism:
#                    baseX=t.solution[radicals[species]]
#                    baseZ=t.Index[0]
#                    baseT=t.solution['T']
#                if t.phi==x and m.replace('forward_','')==t.mechanism:
#                    forX=t.solution[radicals[species]]
#                    forZ=t.Index[0]
#                    forT=t.solution['T']
#            ax1=fig.add_subplot(111)
#            ax1.plot(baseZ,baseX,forZ,forX)
#            ax1.set_ylabel('X')
#            ax2=ax1.twinx()
#            ax2.plot(baseZ,baseT,'r--')
#            ax2.plot(forZ,forT,'r--')
#            ax2.set_ylabel('T',color='r')
#            ax2.set_ylim([2350,2600])
#            plt.title(radicals[species]+' for '+fuels[0]+'_phi='+str(x)+' in '+m.split('\\')[-1].rstrip('.cti').replace('appendedinput','').replace('forward_',''))
#            #plt.savefig(os.getcwd()+'\\figures\\'+fuels[0]+'\\'+m.split('\\')[-1].rstrip('.cti').replace('appendedinput_','').replace('forward_','')+'_phi'+str(x).replace('.','')+radicals[species]+'_radicalProfile'+'.png'
#                       # ,dpi=300,bbox_inches="tight",format='png')

reducedList=[]
for i in mechList:
    reducedList.append(i.rstrip('.cti').split('\\')[-1])
columns=reducedList+['phi','rate']
reactionrateData=pd.DataFrame(columns=columns)
import copy
for i in y:
    finalRates=i.net_rates_of_progress[:,-1]
    
    copyProgress=np.abs(copy.deepcopy(finalRates))
    copyProgress.sort()
    copyProgress=copyProgress.tolist()[-10:]
    indices=[]
    for j in copyProgress:
        indices.append(np.abs(np.array(finalRates)).tolist().index(j))
    reactionLists=[]
    ratevalues=[]
    for j in indices:
        reactionLists.append(i.Index[1][j])
        ratevalues.append(finalRates[j])
    phidata=[i.phi]*10
    tempFrame=pd.DataFrame(columns=[i.mechanism.rstrip('.cti').split('\\')[-1],'phi'])
    tempFrame[i.mechanism.rstrip('.cti').split('\\')[-1]]=reactionLists
    tempFrame['phi']=phidata
    tempFrame['rate']=ratevalues
    reactionrateData=pd.concat([reactionrateData,tempFrame])