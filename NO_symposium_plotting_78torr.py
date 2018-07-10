# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 14:17:50 2017

@author: Mark
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

y=results

linewidths=[2,2,3]
markers=['b-','r--','g:']
fig=plt.figure()
for i in np.arange(len(y)):
    
    
    plt.plot(100*np.array(y[i].Index[0]),y[i].solution['NO']*1000000,markers[i],linewidth=linewidths[i])

title_font = {'fontname':'Arial', 'size':'14'}
plt.annotate(r'$\mathrm{H}_2/$Air, 78 Torr',xy=(1.25,0.05),**title_font)
plt.annotate('Nominal',xy=(1,0.25),xytext=(0.75,0.15),color='b',
             arrowprops=dict(arrowstyle='->',color='b'))

plt.annotate('All Termolecular Rxns',xy=(1,0.31),xytext=(1.15,0.2),color='r',arrowprops=dict(arrowstyle='->',color='r'))
plt.annotate('Only $\mathrm{H+N}_2+\mathrm{O=NH+NO}$',xy=(1.25,0.357),xytext=(1.1,0.4),color='g',arrowprops=dict(arrowstyle='->',color='g'))
plt.xlabel('Axial distance (cm)',**title_font)
plt.ylabel('NO mole fraction, $X_{\mathrm{NO}}$ (ppm)',**title_font)
plt.tick_params(axis='both', which='major', labelsize=12)

datafiles=[os.getcwd()+'\\paperFigs\\harrington_x_78torr.csv']
with open(datafiles[0]) as f:
    profile=f.readlines()
    
for i in np.arange(len(profile)):
    profile[i]=profile[i].rstrip('\n').rstrip(',')+'\n'

with open(datafiles[0],'w') as f:
    f.writelines(profile)
    
data=pd.read_csv(datafiles[0],header=3)
plt.plot(data['z'],data['x'],'kx')
#plt.annotate('Experimental',xy=(0.24,0.255),xytext=(0.0,0.4),arrowprops=dict(arrowstyle='->',color='k'))

#plt.savefig(os.getcwd()+'\\figures\\symposiumFigs\\'+'NO_speciation_78torr.pdf',dpi=1200,bbox_inches='tight')



