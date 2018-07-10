# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 15:09:42 2017

@author: Mark Barbet

"""

import os
os.chdir('C:\Users\HP USER\Google Drive\Burke Group\Codes')
import simulations as sim
import cantera as ct
import numpy as np
import pandas as pd
filename='gri30.xml'
filenames=[os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti',
           os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta001_uscmech.cti',
           os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta001_forward_uscmech.cti']
fuel='C2H4'
phi=0.5
Temps=[1200]
pressures=[29.6077]
resTime=np.linspace(-6,-3,200)
for i in np.arange(len(resTime)):
    resTime[i]=10**resTime[i]
initMolFracs={'C4H10':0.04,'O2':0.26,'HE':0.7}
initMolFracs=[{'C2H4':1.0*phi,'O2':3,'N2':11.28}]

observables=['OH','H2','H2O','CO','CO2']
datafiles=['C:\\Users\\HP USER\\Google Drive\\Burke Group\\Codes\\paperFigs\\ethylene_jsr_restime_oh.csv',
           'C:\\Users\\HP USER\\Google Drive\\Burke Group\\Codes\\paperFigs\\ethylene_jsr_restime_h2.csv',
           'C:\\Users\\HP USER\\Google Drive\\Burke Group\\Codes\\paperFigs\\ethylene_jsr_restime_h2o.csv',
           'C:\\Users\\HP USER\\Google Drive\\Burke Group\\Codes\\paperFigs\\ethylene_jsr_co.csv',
           'C:\\Users\\HP USER\\Google Drive\\Burke Group\\Codes\\paperFigs\\ethylene_jsr_restime_co2.csv']
sens=0

volume=0.00009
results=[]

expConditions=[]
tempresults=[]
for T in Temps:
    for P in pressures:
        for filename in filenames:
            for r in resTime:
                for X in initMolFracs:
                    expConditions.append([T,P,filename,r,X])
for i in expConditions:
    gas=ct.Solution(i[2])
    f={'residenceTime':i[3],'reactorVolume':volume,'pressure':i[1],'conditions':i[4]}

    tempresults.append(sim.multiTemp2(i[2],gas,[i[0]],f,kinetic_sens=0,physical_sens=0,observables=observables,physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=10000))
    gas2=ct.Solution(i[2])

    newCond={}
    for j in gas.species_names:
    #if float(results[0].solution[i])>0.0001:
        newCond[j]=float(tempresults[0].solution[j])
    print(tempresults[0].final_pressure)
    print(i[1])
    #print(newCond)
    #f={'residenceTime':resTime,'reactorVolume':volume,'pressure':1.0,'conditions':newCond} 
    gas2.TPX=i[0],ct.one_atm*i[1],newCond
    results.append(sim.multiTemp2(i[2],gas,[i[0]],f,kinetic_sens=sens,physical_sens=0,observables=observables,physical_params=['T','P'],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=10000,initial_condition_gas=gas2))
    results[-1].add_mechanism(i[2])
    results[-1].assign_phi(phi)
    results[-1].add_residence_time(i[3])
    print(results[-1].final_pressure)
    
mech=[]
tempresults=results
for i in tempresults:
    if i.mechanism not in mech:
        mech.append(i.mechanism)
import matplotlib.pyplot as plt
#print(mech)
for p in np.arange(len(datafiles)):
    with open(datafiles[p]) as f:
        profile=f.readlines()
        
    for i in np.arange(len(profile)):
        profile[i]=profile[i].rstrip('\n').rstrip(',')+'\n'
    
    with open(datafiles[p],'w') as f:
        f.writelines(profile)
    
data=pd.read_csv(datafiles[0],header=3)

count=0
for obs in observables:
    
    data=pd.read_csv(datafiles[count])
    count=count+1
    baseRes=[]
    baseObs=[]
    base=False
    for m in mech:
        tempRes=[]
        tempObs=[]
        
        for i in tempresults:
            #print(i.mechanism)
            if m==i.mechanism and 'appendedinput' not in m:
                baseRes.append(i.residence_time)
                baseObs.append(i.solution[obs])
            elif m==i.mechanism and 'appendedinput' in m:
                base=True
                tempRes.append(i.residence_time)
                tempObs.append(i.solution[obs])
        if base==True:
            plt.figure()
            plt.loglog(baseRes,baseObs,tempRes,tempObs,data['t'],data['x'])
            plt.title('Mole Fraction vs Residence Time for '+obs+':'+' '+m.split('\\')[-1])
            plt.xlabel('Residence Time')
            plt.ylabel('Mole Fraction '+obs)
            plt.savefig(os.getcwd()+'\\figures\\'+fuel+'_'+str(phi)+'_'+m.split('\\')[-1]+'_'+obs+'_molfracs.png',dpi=1200,bbox_inches='tight')
                
                
            
            