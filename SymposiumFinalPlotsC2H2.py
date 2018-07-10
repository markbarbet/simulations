# -*- coding: utf-8 -*-
"""
Created on Sat Nov 04 20:40:54 2017

@author: Mark Barbet
"""

import numpy as np
import pandas as pd
import copy
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

y=[]
for i in results:
    if i!='failed':
        y.append(i)

finalsubset=[]
finalmechs=['FFCM1.cti','appendedinput_beta001_forward_FFCM1.cti','appendedinput_beta001_forward_top30_FFCM1.cti']

for i in y:
    if i.mechanism.split('\\')[-1] in finalmechs:
        finalsubset.append(i)

finalFlamespeeds=pd.DataFrame(columns=finalmechs+['phi'])

for i in np.arange(len(finalmechs)):
    
    tempflames=[]
    tempphi=[]
    for j in finalsubset:
        #print(j.mechanism.split('\\')[-1])
        if j.mechanism.split('\\')[-1]==finalmechs[i]:
            tempflames.append(j.solution['u'][0]*100)
            tempphi.append(j.phi)
    finalFlamespeeds[finalmechs[i]]=tempflames
    finalFlamespeeds['phi']=tempphi

deviations=pd.DataFrame()
deviations['phi']=tempphi
for i in finalmechs:
    if 'appendedinput' not in i:
        base=i
    elif 'appendedinput' in i:
        deviations[i]=np.divide(finalFlamespeeds[i]-finalFlamespeeds[base],finalFlamespeeds[base])
#from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})
markerlist=['b-','r--','g:']
fig=plt.figure(figsize=(4,4))
gs = gridspec.GridSpec(2, 1,height_ratios=[3,1],wspace=0.025,hspace=0.1)
ax1=plt.subplot(gs[0])
ax2=plt.subplot(gs[1])
#f,(ax1,ax2)=plt.subplots(1,2,sharex=False,sharey=False,figsize=(15,15))
for i in np.arange(len(finalmechs)):
    
    ax1.plot(finalFlamespeeds['phi'],finalFlamespeeds[finalmechs[i]],markerlist[i])
    if 'appendedinput' in finalmechs[i]:
        ax2.plot(deviations['phi'],100*deviations[finalmechs[i]],markerlist[i])
axis_font = {'fontname':'Arial', 'size':'14'}
ax1.set_ylabel('Flame speed, $S_u^0$ (cm/s)',**axis_font)
ax1.text(0.57,107,'$\mathrm{C}_2\mathrm{H}_2$',**axis_font)
#ax1.text(3.75,130,b'Nominal Model')
ax1.annotate('Nominal',xy=(1.6,90),xytext=(1.4,105),color='b',arrowprops=dict(arrowstyle='->',color='b'))
ax1.annotate('All Termolecular Rxns',xy=(0.7,65),xytext=(0.8,60),color='r',arrowprops=dict(arrowstyle='->',color='r'))
ax1.annotate('Only Influential Rxns',xy=(1.2,108),xytext=(0.9,70),color='g',
             arrowprops=dict(arrowstyle='->',color='g'))
ax1.tick_params(axis='both',direction='in')
ax1.set_xticks([0.75,1.0,1.25,1.5,1.75])
ax1.set_xticklabels('')
#ax2.yaxis.tick_right()
ax2.tick_params(axis='both',direction='in')
ax2.yaxis.set_label_position('right')
#ax2.set_ylabel('Relative Difference',**axis_font)
ax2.annotate('% Deviation from Nominal',xy=(1.05,-8.5))
#ax1.set_xlabel('Equivalence Ratio, $\phi$')
#ax2.set_xlabel('Equivalence Ratio, $\phi$')
fig.text(0.5, 0.04, 'Equivalence ratio, $\phi$', ha='center',**axis_font)
fig.subplots_adjust(left=None, bottom=0.15, right=None, top=None, wspace=None, hspace=None)

#datafiles=[os.getcwd()+'\\paperFigs\\c2h2_flamespeeds.csv']
#with open(datafiles[0]) as f:
#    profile=f.readlines()
#    
#for i in np.arange(len(profile)):
#    profile[i]=profile[i].rstrip('\n').rstrip(',')+'\n'
#
#with open(datafiles[0],'w') as f:
#    f.writelines(profile)
#    
#data=pd.read_csv(datafiles[0],header=3)
#ax1.plot(data['phi'],data['f'],'kx')
#ax1.annotate('Experiment',xy=(data['phi'][6],data['f'][6]),xytext=(4.5,120),arrowprops=dict(arrowstyle='->',color='k'))

#plt.savefig(os.getcwd()+'\\figures\\symposiumFigs\\'+'C2H2_flamespeed_fix.pdf',dpi=1200,bbox_inches='tight')
