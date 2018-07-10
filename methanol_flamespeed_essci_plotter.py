# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 18:30:14 2018

@author: Mark Barbet
"""
#Code to analyze methanol from the dryer mechanism for methanol flame speeds




import numpy as np
import pandas as pd
import copy
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#y=[]
#for i in results:
#    if i!='failed':
#        y.append(i)
#print('Hi!')

#Section to construct the final datasets for plotting
nominal=[]
modified=[]
phin=[]
phim=[]
for i in results:
    if 'efficiency_manipulated:True' in i.tags:
        modified.append(i.solution['u'][0])
        phim.append(i.phi)
    elif 'efficiency_manipulated:False' in i.tags:
        nominal.append(i.solution['u'][0])
        phin.append(i.phi)
difference=[]
for j in np.arange(len(phin)):
    if phin[j] in phim:
        index=phim.index(phin[j])
        difference.append(100.0*np.divide(modified[index]-nominal[j],nominal[j]))
axis_font = {'fontname':'Arial', 'size':'13'}         
axis_font2= {'fontname':'Arial','size':'11'}
markerlist=['b-','r--','g:']
fig=plt.figure(figsize=(4,4))
gs = gridspec.GridSpec(2, 1,height_ratios=[3,1],wspace=0.025,hspace=0.1)
ax1=plt.subplot(gs[0])
ax2=plt.subplot(gs[1])
#f,(ax1,ax2)=plt.subplots(1,2,sharex=False,sharey=False,figsize=(15,15))
ax1.plot(phin,nominal,markerlist[0])
ax1.plot(phim,modified,markerlist[1])
#ax2.yaxis.tick_right()
ax2.plot(phim,difference,markerlist[1])
ax1.set_ylabel('Flame speed, $S_u^0$ (cm/s)',**axis_font)
ax2.yaxis.set_label_position("left")
#ax2.set_ylabel('Relative Percent Difference',**axis_font)
ax1.set_xticks([0.7,1.0,1.3,1.6,1.9,2.2])
ax1.set_xticklabels('')

ax2.set_xticks([0.7,1.0,1.3,1.6,1.9,2.2])
ax2.set_yticks([0,20.,40.])
ax1.annotate('      $\mathrm{CH_3OH}$/air flames\n$p=1$ atm, $T=298$ K',xy=(1.5,0.44),color='k',**axis_font2)
ax2.set_xlabel('Equivalence Ratio, $\phi$',**axis_font)
ax1.annotate('Nominal Model',xy=(1.9,0.15),xytext=(0.7,0.1),color='b',arrowprops=dict(arrowstyle='->',color='b'))
ax1.annotate('Modified Model',xy=(2.2,0.12),xytext=(1.6,0.35),color='r',arrowprops=dict(arrowstyle='->',color='r'))
ax1.tick_params(axis='both',direction='in',labelsize=12)
ax2.tick_params(axis='both',direction='in',labelsize=12)
ax2.annotate('% Deviation from Nominal',xy=(0.60,35),**axis_font)
#plt.savefig(os.getcwd()+'\\figures\\essci2018\\methanol_flamespeed.pdf',dpi=1200,bbox_inches='tight')
#for i in np.arange(len(finalmechs)):
#    
#    ax1.plot(finalFlamespeeds['phi'],finalFlamespeeds[finalmechs[i]],markerlist[i])
#    if 'appendedinput' in finalmechs[i]:
#        ax2.plot(deviations['phi'],100*deviations[finalmechs[i]],markerlist[i])
#axis_font = {'fontname':'Arial', 'size':'14'}
#ax1.set_ylabel('Flame speed, $S_u^0$ (cm/s)',**axis_font)
#ax1.text(0.57,107,'$\mathrm{C}_2\mathrm{H}_2$',**axis_font)
##ax1.text(3.75,130,b'Nominal Model')
#ax1.annotate('Nominal',xy=(1.6,90),xytext=(1.4,105),color='b',arrowprops=dict(arrowstyle='->',color='b'))
#ax1.annotate('All Termolecular Rxns',xy=(0.7,65),xytext=(0.8,60),color='r',arrowprops=dict(arrowstyle='->',color='r'))
#ax1.annotate('Only Influential Rxns',xy=(1.2,108),xytext=(0.9,70),color='g',
#             arrowprops=dict(arrowstyle='->',color='g'))
#ax1.tick_params(axis='both',direction='in')
#ax1.set_xticks([0.75,1.0,1.25,1.5,1.75])
#ax1.set_xticklabels('')
##ax2.yaxis.tick_right()
#ax2.tick_params(axis='both',direction='in')
#ax2.yaxis.set_label_position('right')
##ax2.set_ylabel('Relative Difference',**axis_font)
#ax2.annotate('% Deviation from Nominal',xy=(1.05,-8.5))
##ax1.set_xlabel('Equivalence Ratio, $\phi$')
##ax2.set_xlabel('Equivalence Ratio, $\phi$')
#fig.text(0.5, 0.04, 'Equivalence ratio, $\phi$', ha='center',**axis_font)
#fig.subplots_adjust(left=None, bottom=0.15, right=None, top=None, wspace=None, hspace=None)

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
