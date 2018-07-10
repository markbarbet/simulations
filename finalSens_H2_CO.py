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

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <=rel_tol
mechs=['appendedinput_beta001_forward_FFCM1.cti']
philist=[0.6,1.0,5.0]
num_oldrxns=0
width=0.2
step=0
for i in results:
    if i.mechanism.split('\\')[-1]=='FFCM1.cti':
        num_oldrxns=len(i.Index[1])
fig,ax=plt.subplots()


subresults=[]
for k in results:
    for j in np.arange(len(philist)):
        if k.mechanism.split('\\')[-1] in mechs and isclose(philist[j],k.phi):
            #print('hi')
            subresults.append(k)
colors=['b','r','g']   

fsens=subresults[0].flamespeed_sens[num_oldrxns:]['Su']
rxns=subresults[0].Index[1][num_oldrxns:]
absVal=np.abs(fsens)
temp=pd.DataFrame(columns=['sens','abs','rxn'])
temp['sens']=fsens
temp['abs']=absVal
temp['rxn']=rxns
sortedvals=temp.sort_values(by=['abs'],ascending=False)
topsens=sortedvals['sens'][0:10]
toprxns=sortedvals['rxn'][0:10]
#indexing=np.arange(len(toprxns))
        
#ax.barh(indexing+step,topsens,width,color=colors[0])
#step=step+0.2

#ax.set_yticks(np.arange(len(toprxns)))
#ax.set_yticklabels(toprxns)

fsens2=subresults[1].flamespeed_sens[num_oldrxns:]['Su']
rxns2=subresults[1].Index[1][num_oldrxns:]
absVal2=np.abs(fsens2)
temp2=pd.DataFrame(columns=['sens','abs','rxn'])
temp2['sens']=fsens2
temp2['abs']=absVal2
temp2['rxn']=rxns2
sortedvals2=temp2.sort_values(by=['abs'],ascending=False)
topsens2=sortedvals2['sens'][0:10]
toprxns2=sortedvals2['rxn'][0:10]

fsens3=subresults[2].flamespeed_sens[num_oldrxns:]['Su']
rxns3=subresults[2].Index[1][num_oldrxns:]
absVal3=np.abs(fsens3)
temp3=pd.DataFrame(columns=['sens','abs','rxn'])
temp3['sens']=fsens3
temp3['abs']=absVal3
temp3['rxn']=rxns3
sortedvals3=temp3.sort_values(by=['abs'],ascending=False)
topsens3=sortedvals3['sens'][0:10]
toprxns3=sortedvals3['rxn'][0:10]

    