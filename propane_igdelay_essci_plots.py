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
plt.annotate('Modified Model',xy=(0.71,10**-3.4),xytext=(0.685,10**-3.2),color='r',arrowprops=dict(arrowstyle='->',color='r'))
plt.annotate('Nominal Model',xy=(0.71,10**-3.44),xytext=(0.73,10**-3.6),color='b',arrowprops=dict(arrowstyle='->',color='b'))
plt.annotate(r'$\mathrm{C_3H_8}$/$\mathrm{O_2}$/$\mathrm{Ar}$ = 1.6/8/90.4',xy=(0.685,10**-2.87),color='k',**axis_font)
plt.annotate('Shock tube ignition delay times',xy=(0.685,10**-2.78),color='k',**axis_font)
plt.tick_params(axis='both',direction='in',labelsize=12)
plt.xticks([0.69,0.71,0.73,0.75,0.77])
plt.savefig(os.getcwd()+'\\figures\\essci2018\\propane_igdelay.pdf',dpi=1200,bbox_inches='tight')
plt.savefig(os.getcwd()+'\\figures\\essci2018\\propane_igdelay.png',dpi=1200,bbox_inches='tight')