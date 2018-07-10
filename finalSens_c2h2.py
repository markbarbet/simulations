# -*- coding: utf-8 -*-
"""
Created on Mon Nov 06 10:19:27 2017

@author: Skoron
"""

import numpy as np
import pandas as pd
import os
import copy
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <=rel_tol
mechs=['appendedinput_beta001_forward_FFCM1.cti']
philist=[0.6,1.0,1.8]
philist=[1.8,1.0,0.6]
philist=[0.6,1.0,1.8]
num_oldrxns=0
width=0.3
step=0
for i in results:
    if i.mechanism.split('\\')[-1]=='FFCM1.cti':
        num_oldrxns=len(i.Index[1])
fig,ax=plt.subplots(figsize=(7,5))

axis_font = {'fontname':'Arial', 'size':'14'}    
subresults=[]

for j in np.arange(len(philist)):
    for k in results:
        if k.mechanism.split('\\')[-1] in mechs and isclose(philist[j],k.phi):
            print(j)
            subresults.append(k)
phicond=['$\phi=0.6$','$\phi=1.0$','$\phi=1.8$']
#phicond=phicond[::-1]
colors=['b','r','w']   
fsens=subresults[0].flamespeed_sens[num_oldrxns:]['Su']
rxns=subresults[0].Index[1][num_oldrxns:]
absVal=np.abs(fsens)
temp=pd.DataFrame(columns=['sens','abs','rxn'])
temp['sens']=fsens
temp['abs']=absVal
temp['rxn']=rxns
sortedvals=temp.sort_values(by=['abs'],ascending=False)
topsens=sortedvals['sens'][0:10].tolist()
toprxns=sortedvals['rxn'][0:10].tolist()
indexing=np.arange(len(toprxns))
ax.barh(indexing-step+0.2,topsens[::-1],width,color=colors[0],label=phicond[0],edgecolor='k', linewidth=1)
labels=[]
for l in toprxns:
    labels.append(l.split('=')[0].rstrip('<')+'='+l.split('=')[1].lstrip('>'))
maxlen=0
#for l in np.arange(len(labels)):
#    if len(labels[l])>maxlen:
#        
#        maxlen=len(labels[l])
#for l in np.arange(len(labels)):
#    if len(labels[l])<maxlen:
#        padding=maxlen-len(labels[l])
#        labels[l]=labels[l]+' '*padding

#labels=['CO+2H=$\mathrm{H}_2$','2H+$\mathrm{O}_2(+\mathrm{M})=\mathrm{H}_2+\mathrm{O}_2(+\mathrm{M})$',]
ax.set_yticks(np.arange(len(toprxns)))
#ax.set_yticklabels(labels[::-1],position=(0.05,0),horizontalalignment='left')
ax.set_yticklabels(''*len(toprxns),position=(0.05,0),horizontalalignment='left')
equations=['$\mathrm{H+O}_2$+H(+M)=2OH(+M)','H+$\mathrm{C}_2\mathrm{H}_2$+$\mathrm{O}_2$(+M)=$\mathrm{CH}_2$O+HCO(+M)','H+O+$\mathrm{O}_2$(+M)=$\mathrm{O}_2$+OH(+M)',
           'H+$\mathrm{O}_2$+H(+M)=$\mathrm{H}_2$+$\mathrm{O}_2$(+M)','H+$\mathrm{O}_2$+OH(+M)=$\mathrm{H}_2$O+$\mathrm{O}_2$(+M)','H+$\mathrm{C}_2\mathrm{H}_2$+$\mathrm{O}_2$(+M)=$\mathrm{CH}_2$CHO+O(+M)',
           'H+$\mathrm{O}_2$+OH(+M)=$\mathrm{H}_2$O+$\mathrm{O}_2$(+M)','H+$\mathrm{C}_2\mathrm{H}_2$+H(+M)=$\mathrm{H}_2$+$\mathrm{H}_2$CC(+M)',
           'H+$\mathrm{C}_2\mathrm{H}_2$+OH(+M)=$\mathrm{C}_2\mathrm{H}_2$+$\mathrm{H}_2$O(+M)','H+CO+OH=CO+$\mathrm{H}_2$O']
xys=[(0.0005,-0.25),(0.0005,1-0.25),(0.0005,2-0.25),(0.0005,3-0.25),(-0.018,4-0.25),(0.0005,5-0.25),(0.0005,6-0.25),(0.0005,7-0.25),(0.0005,8-0.25),(0.0005,9-0.25)]
xys=xys[::-1]
ax.set_xlim(right=0.018)
for i in np.arange(len(xys)):
    ax.annotate(equations[i],xy=xys[i])
hatches=['xxxx','/////']
for i in np.arange(len(philist)):
    if i>0:
        step=step+width
        fsens1=subresults[i].flamespeed_sens[num_oldrxns:]['Su']
        rxns1=subresults[i].Index[1][num_oldrxns:]
        absVal1=np.abs(fsens1)
        temp1=pd.DataFrame(columns=['sens','abs','rxn'])
        temp1['sens']=fsens1
        temp1['abs']=absVal1
        temp1['rxn']=rxns1
        sortedvals1=temp1.sort_values(by=['abs'],ascending=False)
        
        topsens1=sortedvals1['sens'][0:10]
        toprxns1=sortedvals1['rxn'][0:10]
        indexing1=np.arange(len(toprxns1))
        senstoplot=sortedvals1[sortedvals1['rxn'].isin(toprxns)]
        possiblesens=senstoplot['sens'].tolist()
        possiblerxns=senstoplot['rxn'].tolist()
        matched=[False]*10
        matched2=[False]*len(senstoplot['rxn'])
        finalsensitivities=[]
        for j in np.arange(len(toprxns)):
            for k in np.arange(len(senstoplot['rxn'])):
                if toprxns[j]==senstoplot['rxn'][k] and matched[j]==False and matched2[k]==False:
                    matched[j]=True
                    matched2[k]=True
                    finalsensitivities.append(senstoplot['sens'][k])
        
        ax.barh(indexing-step+0.2,finalsensitivities[::-1],width,color=colors[i],label=phicond[i],hatch=hatches[i-1],edgecolor='k', linewidth=1)
    
ax.tick_params(axis='y',direction='in')
ax.tick_params(axis='x',labelsize=12,direction='in')
ax.annotate('$\mathrm{C}_2\mathrm{H}_2$',xy=(-0.018,-0.5),**axis_font)
ax.set_xlabel(r'Sensitivity coefficient, $\frac{\partial\mathrm{ln}S_u^0}{\partial\mathrm{ln}k}$',**axis_font)
handles,labels = ax.get_legend_handles_labels()
handles=[handles[0],handles[1],handles[2]]
labels=[labels[0],labels[1],labels[2]]

ax.legend(handles,labels,loc=3,bbox_to_anchor=(0.01,0.1))
#ax.legend(phicond,loc=3,bbox_to_anchor=(0.01,0.1))
#ax.add_line(Line2D([0.5, 0.5], [0, 1], transform=ax.transAxes,
#                  linewidth=1, color='k'))
ax.axvline(x=0, color='k',linewidth=1.0)
plt.savefig(os.getcwd()+'\\figures\\symposiumFigs\\'+'C2H2_flamespeed_sens_symposium_fix_May2018.pdf',dpi=1200,bbox_inches='tight')
#    ax.barh(indexing+step,topsens,width,color=colors[0])
#    step=step+0.2
#    
#    ax.set_yticks(np.arange(len(toprxns)))
#    ax.set_yticklabels(toprxns)

#fsens2=subresults[1].flamespeed_sens[num_oldrxns:]['Su']
#rxns2=subresults[1].Index[1][num_oldrxns:]
#absVal2=np.abs(fsens2)
#temp2=pd.DataFrame(columns=['sens','abs','rxn'])
#temp2['sens']=fsens2
#temp2['abs']=absVal2
#temp2['rxn']=rxns2
#sortedvals2=temp2.sort_values(by=['abs'],ascending=False)
#topsens2=sortedvals2['sens'][0:10]
#toprxns2=sortedvals2['rxn'][0:10]
#
#fsens3=subresults[2].flamespeed_sens[num_oldrxns:]['Su']
#rxns3=subresults[2].Index[1][num_oldrxns:]
#absVal3=np.abs(fsens3)
#temp3=pd.DataFrame(columns=['sens','abs','rxn'])
#temp3['sens']=fsens3
#temp3['abs']=absVal3
#temp3['rxn']=rxns3
#sortedvals3=temp3.sort_values(by=['abs'],ascending=False)
#topsens3=sortedvals3['sens'][0:10]
#toprxns3=sortedvals3['rxn'][0:10]

    