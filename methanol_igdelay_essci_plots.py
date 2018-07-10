# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 17:01:07 2018

@author: Mark Barbet
"""
#File to plot graphs for methanol ignition delays

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

datafiles=[os.getcwd()+'\\paperFigs\\methanol_ig_delay.csv']
with open(datafiles[0]) as f:
    profile=f.readlines()
    
for i in np.arange(len(profile)):
    profile[i]=profile[i].rstrip('\n').rstrip(',')+'\n'

with open(datafiles[0],'w') as f:
    f.writelines(profile)
    
data=pd.read_csv(datafiles[0],header=3)

plt.figure(figsize=(4,4))
colors=['b-','r--','k','g','c','m']
for i in np.arange(len(mechanism)):
    plt.semilogy(1000.0/T,delays[i],colors[i])
    
#plt.semilogy(data['x'],10**data['logt'])
axis_font = {'fontname':'Arial', 'size':'13'}         

plt.xlabel('1000 K/T',**axis_font)
plt.ylabel(r'$\tau_{\mathrm{IDT}}$ / s',**axis_font)
plt.annotate('Modified Model',xy=(0.9,10**-3.60),xytext=(0.9,10**-4.0),color='r',arrowprops=dict(arrowstyle='->',color='r'))
plt.annotate('Nominal Model',xy=(0.922,10**-3.25),xytext=(0.85,10**-3.125),color='b',arrowprops=dict(arrowstyle='->',color='b'))
plt.annotate(r'$\mathrm{CH_3OH}$/$\mathrm{O_2}$/$\mathrm{N_2}$ = 5/15/80',xy=(0.805,10**-2.85),color='k',**axis_font)
plt.annotate('Shock tube ignition delay times',xy=(0.805,10**-2.7),color='k',**axis_font)
plt.tick_params(axis='both',direction='in',labelsize=12)
plt.xticks([0.8,0.85,0.9,0.95,1.0])
plt.savefig(os.getcwd()+'\\figures\\essci2018\\methanol_igdelay.pdf',dpi=1200,bbox_inches='tight')
plt.savefig(os.getcwd()+'\\figures\\essci2018\\methanol_igdelay.png',dpi=1200,bbox_inches='tight')