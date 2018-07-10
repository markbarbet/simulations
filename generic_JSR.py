# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 18:25:42 2017

@author: Mark Barbet
"""

#JSR loop to do sensitivity analysis on steady state mol fracs


import cantera as ct
import numpy as np
import os
import simulations as sim


os.chdir('C:\Users\HP USER\Google Drive\Burke Group\Codes')

#########################################################################################
#########################################################################################
#Change inputs here

Temps=np.linspace(800,1100,4)
filenames=[os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1.cti','gri30.xml']

pressures=[1]
phi=1
fuel='CH4'
resTime=[2]
conditions=[{fuel:1,'O2':1,'N2':3.76}]
observables=['O2']
volume=0.00009

sens=1   #Edit to 1 to run sensitivities

#End inputs here
#########################################################################################
#########################################################################################


expConditions=[]

for filename in filenames:
    for T in Temps:
        for P in pressures:
            for r in resTime:
                for X in conditions:
                    expConditions.append([T,P,filename,r,X])
tempresults=[]
results=[]
for i in expConditions:
    gas=ct.Solution(i[2])
    f={'residenceTime':i[3],'reactorVolume':volume,'pressure':i[1],'conditions':i[4]}
    if sens==1:
        tempresults.append(sim.multiTemp2(i[2],gas,[i[0]],f,kinetic_sens=0,physical_sens=0,observables=observables,physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=10000))
        gas2=ct.Solution(i[2])
        newCond=0
        newCond={}
        for j in gas.species_names:
            #if float(results[0].solution[i])>0.0001:
            newCond[j]=float(tempresults[-1].solution[j])
        f={'residenceTime':i[3],'reactorVolume':volume,'pressure':1.0,'conditions':newCond}
        gas2.TPX=i[0],ct.one_atm*i[1],newCond
        results.append(sim.multiTemp2(i[2],gas,[i[0]],f,kinetic_sens=sens,physical_sens=0,observables=observables,physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=10000,initial_condition_gas=gas2))
#        print(tempresults[0].final_pressure)
#        print(i[1])
    #print(newCond)
    #f={'residenceTime':resTime,'reactorVolume':volume,'pressure':1.0,'conditions':newCond} 
    else:
        #gas2.TPX=i[0],ct.one_atm*i[1],newCond
        results.append(sim.multiTemp2(i[2],gas,[i[0]],f,kinetic_sens=sens,physical_sens=0,observables=observables,physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=10000))
    results[-1].add_mechanism(i[2])
    results[-1].add_fuel(fuel)
    results[-1].assign_phi(phi)
    results[-1].add_residence_time(i[3])
    print(results[-1].final_pressure)