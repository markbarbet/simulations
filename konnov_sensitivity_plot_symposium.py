# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 15:40:09 2017

@author: Mark
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

p=results

finalsubset=[]
for i in p:
    if i.mechanism.split('\\')[-1] in ['appendedinput_beta1_forward_chem_onlyNO_H2.cti','appendedinput_beta1_forward_chem_onlyNO.cti']:
        finalsubset.append(i)
axis_font = {'fontname':'Arial', 'size':'13'}         
fig,ax=plt.subplots(figsize=(6,2))

pressures=['38 Torr','78 Torr']

dframe=pd.DataFrame(columns=['rxn']+['sens']+['abs'])
dframe['rxn']=finalsubset[1].Index[1]
dframe['sens']=finalsubset[1].k_sens[-1,:,:]
dframe['abs']=np.abs(dframe['sens'])
sortedvals=dframe.sort_values(by=['abs'],ascending=False)
finalvals=sortedvals[0:3]
sens=finalvals['sens'].tolist()
rxns=finalvals['rxn'].tolist()

step=0.45
width=0.45

dframe2=pd.DataFrame(columns=['rxn']+['sens']+['abs'])
dframe2['rxn']=finalsubset[0].Index[1]
dframe2['sens']=finalsubset[0].k_sens[-1,:,:]
dframe2['abs']=np.abs(dframe2['sens'])
sortedvals2=dframe2.sort_values(by=['abs'],ascending=False)
#finalvals2=sortedvals2[0:15]
matched=[False]*10
matched2=[False]*len(sortedvals2)
possiblesens=sortedvals2['sens'].tolist()
possiblerxns=sortedvals2['rxn'].tolist()
sens2=[]
rxns2=[]
for i in np.arange(len(sortedvals[0:3])):
    #print(i)
    for j in np.arange(len(sortedvals2)):
        #print(j)
        if matched[i]==False and matched2[j]==False and rxns[i]==possiblerxns[j]:
            #print(rxns[i],sortedvals2['rxn'][j],i)
            sens2.append(possiblesens[j])
            rxns2.append(possiblerxns[j])
            matched[i]=True
            matched2[j]=True
labels=[]
for i in rxns:
    labels.append(i.split('=')[0].rstrip('<')+'='+i.split('=')[1].lstrip('>'))
labels=['H+O$_2$(+M) = HO$_2$(+M)','H+N$_2$+O=NH+NO','H$_2$+O=H+OH']
labels=labels[::-1]
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(['']*len(labels))
indexing=np.arange(len(rxns))
ax.barh(indexing-0.25,sens[::-1],width,color='b',label=pressures[1],edgecolor='k', linewidth=1)
ax.barh(indexing+step-0.25,sens2[::-1],width,color='w',label=pressures[0],hatch='////',edgecolor='k', linewidth=1)

handles,labels1 = ax.get_legend_handles_labels()
handles=[handles[1],handles[0]]
labels1=[labels1[1],labels1[0]]

#ax.legend(pressures[::-1],loc=1,bbox_to_anchor=(0.74,0.37))
ax.legend(handles,labels1,loc=1,bbox_to_anchor=(0.74,0.37))
ax.set_xlabel(r'Sensitivity coefficient, $\frac{\partial\mathrm{ln}[NO]}{\partial\mathrm{ln}k}$',**axis_font)
ax.axvline(x=0, color='k',linewidth=1.0)
ax.annotate('$\mathrm{H}_2/$Air:$\phi=1.5$',xy=(0.36,-0.25),**axis_font)
xys=[(-0.7,0.1),(-0.7,0.9),(0.1,1.9),(-0.7,2.9),(0.1,3.9)]
xys=xys[0:3]
for i in np.arange(len(xys)):
    ax.annotate(labels[i],xy=xys[i])
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both',direction='in')
plt.savefig(os.getcwd()+'\\figures\\symposiumFigs\\'+'konnov_harrington_sens.pdf',dpi=1200,bbox_inches='tight')