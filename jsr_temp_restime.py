# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 13:57:24 2017

@author: HP USER
"""

'''
Some code to attempt to generate the weird Tianfeng Liu plots
'''

import cantera as ct
import numpy as np
import os
import simulations as sim


os.chdir('C:\\Users\\HP USER\\Google Drive\\Burke Group\\codes')
filenames=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
           os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_FFCM1.cti',
           os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1.cti']
#filenames=['gri30.xml']
fuel='C2H6'
phi=1.0
conditions=[{fuel:0.1*phi,'O2':0.35,'AR':0.55}]
sens=1
T_inlet=[1000]
pressures=[0.5]
volume=0.00009
resTimes=np.linspace(0,-6,100)
#resTimes=np.array([-4.0,-3.0,-1.0,0.0])
for i in np.arange(len(resTimes)):
    resTimes[i]=10**resTimes[i]
sensOption=np.zeros(len(resTimes))
for i in np.arange(len(sensOption)):
    if i%20==0:
        sensOption[i]=1
initialCond=[{fuel:0.1*phi,'O2':0.35,'H':0.55}]
initialRadicals=[True,False]
expConditions=[]
for T in T_inlet:
    for P in pressures:
        for r in resTimes:
            for X in conditions:
                for f in filenames:
                    for b in initialRadicals:
                        expConditions.append([T,P,r,X,f,b])
results=[]  
tempresults=[]                           
for f in filenames:
                 
    newCond={}
    initialGas=ct.Solution(f)
    gas=ct.Solution(filenames[0])
    initialGas.TPX=T_inlet[0],pressures[0],conditions[0]
    gas.TPX=T_inlet[0],pressures[0],conditions[0]
    c={'residenceTime':resTimes[0],'reactorVolume':volume,'pressure':pressures[0],'conditions':conditions[0]}
    results.append(sim.multiTemp2(f,gas,T_inlet,c,kinetic_sens=sens,physical_sens=0,physical_params=['T','P'],energycon='on',pressureValveCoefficient=1,maxsimulationTime=10000))
    results[-1].add_residence_time(resTimes[0])
    results[-1].add_mechanism(f)
    results[-1].add_fuel(fuel)
    results[-1].assign_phi(phi)
    
    sensConditions={}
    dk=0.01
    tempsens=np.zeros(gas.n_reactions)
    for k in gas.species_names:
        sensConditions[k]=float(results[-1].solution[k])
    sensGas=ct.Solution(filenames[0])
    sensGas.TPX=T_inlet[0],pressures[0],conditions[0]
    sensResult=[]
    for k in np.arange(gas.n_reactions):
        
        sensGas.set_multiplier(1.0)
        sensGas.set_multiplier(1+dk,k)
        senstempSoln=sim.multiTemp2(f,sensGas,T_inlet,c,kinetic_sens=0,physical_sens=0,physical_params=['T','P'],energycon='on',pressureValveCoefficient=1,maxsimulationTime=10000)
        sensResult.append(senstempSoln.solution['temperature'])
        sensResult[-1]=np.log10(sensResult[-1])-np.log10(results[-1].solution['temperature'])
        sensResult[-1]=sensResult[-1]/dk
    results[-1].add_temperature_sensitivity(sensResult)
    for i in np.arange(len(resTimes)):
        if i>0:
            for j in gas.species_names:
                newCond[j]=float(results[-1].solution[j])
            gas2=ct.Solution(f)
            gas2.TPX=results[-1].solution['temperature'],pressures[0]*ct.one_atm,newCond
            
            c={'residenceTime':resTimes[i],'reactorVolume':volume,'pressure':pressures[0],'conditions':conditions[0]}
            results.append(sim.multiTemp2(f,gas,T_inlet,c,kinetic_sens=sensOption[i],physical_sens=0,physical_params=['T','P'],energycon='on',pressureValveCoefficient=1,maxsimulationTime=10000,initial_condition_gas=gas2))
            results[-1].add_residence_time(resTimes[i])
            results[-1].add_mechanism(f)
            results[-1].add_fuel(fuel)
            results[-1].assign_phi(phi)
            if i%20==0:
                sensConditions={}
                dk=0.01
                tempsens=np.zeros(gas.n_reactions)
                for k in gas.species_names:
                    sensConditions[k]=float(results[-1].solution[k])
                sensGas=ct.Solution(filenames[0])
                sensGas.TPX=T_inlet[0],pressures[0],conditions[0]
                sensResult=[]
                for k in np.arange(gas.n_reactions):
            
                    sensGas.set_multiplier(1.0)
                    sensGas.set_multiplier(1+dk,k)
                    senstempSoln=sim.multiTemp2(f,sensGas,T_inlet,c,kinetic_sens=0,physical_sens=0,physical_params=['T','P'],energycon='on',pressureValveCoefficient=1,maxsimulationTime=10000)
                    sensResult.append(senstempSoln.solution['temperature'])
                    sensResult[-1]=np.log10(sensResult[-1])-np.log10(results[-1].solution['temperature'])
                    sensResult[-1]=sensResult[-1]/dk
                results[-1].add_temperature_sensitivity(sensResult)
r=[]
T=[]
for i in results:
    r.append(i.residence_time)
    T.append(i.solution['temperature'])
    
import matplotlib.pyplot as plt
plt.semilogx(r,T)
#runfile('C:/Users/HP USER/Google Drive/Burke Group/Codes/jsr_restime_postprocess.py', wdir='C:/Users/HP USER/Google Drive/Burke Group/Codes')