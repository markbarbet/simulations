# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 16:58:04 2017

@author: HP USER
"""
import matplotlib.pyplot as plt
import cantera as ct
import numpy as np
import JSR_funcs as jsr



gas = ct.Solution('C:\\Users\\HP USER\\Desktop\\kevin\\nbutane\\butane.cti')
speciesName = np.array(['C4H10','O2','B2CO'])
Temps = np.linspace(550,800,50)
component_final_X=np.zeros((len(Temps),len(speciesName)),float)




x,y=jsr.JSR_isothermal_stdst(Temps,gas,ct.one_atm,{'C4H10': 4, 'O2': 26, 'HE': 70},6,30.5*(1e-2)**3)


state=np.zeros(len(Temps))
for k in np.arange(len(Temps)):
    state = np.hstack([x[k].thermo.P, x[k].mass, 
                   x[k].volume, x[k].T, x[k].thermo.X])
    for j in np.arange(len(speciesName)):
        component_final_X[k,j] = state[x[k].component_index(speciesName[j])+1]
        

plt.clf()
a=plt.figure(1)    
plt.plot(Temps,component_final_X[:,0])
plt.xlabel('Temperature (K)')
plt.ylabel(r'Mole Fraction: '+'Butane')