# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 15:20:28 2017

@author: Mark Barbet
"""

#Burner Flame code for searching for interesting termolecular reactions in stabilized flames

import numpy as np
import cantera as ct
import pandas as pd
import os
import simulations as sim
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
import matplotlib.pyplot as plt


phi=[0.7,1,1.3]
phi=np.arange(0.6,1.3,0.2)
phi=[]
#phi=[1]
conditions=[{'C4H6':0.0786,'O2':0.231,'AR':0.2393}]
temps=[566]
pressures=[0.0395]
mdots=[0.0439]
datafiles=[os.getcwd()+'\\paperFigs\\18_T_data.txt.csv']
mechanisms=[os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti',
            os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta1_forward_uscmech.cti',
            os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta1_uscmech.cti']

global conditionsTups
width=0.03
conditionsTups=[]




for i in datafiles:
    with open(i) as f:
        profile=f.readlines()
    
    for j in np.arange(len(profile)):
        profile[j]=profile[j].rstrip('\n').rstrip(',')+'\n'
    
    with open(i,'w') as f:
        f.writelines(profile)
data=[]
for i in datafiles:
    data.append(pd.read_csv(i,header=3))
    
for i in np.arange(len(conditions)):
    for j in mechanisms:
        conditionsTups.append([conditions[i],j,mdots[i],data[i],pressures[i],temps[i]])

def solver(conditionsTup):
    try:
        #print(conditionsTup[0])
        gas=ct.Solution(conditionsTup[1])
        grid=np.linspace(0,0.03,10)
        mdot=conditionsTup[2]
        gas.TPX=conditionsTup[5],conditionsTup[4]*ct.one_atm,conditionsTup[0]
        print(conditionsTup[0])
        results=sim.burner_flame(gas,grid,mdot,data=conditionsTup[3],kinetic_sens=1)
        results.add_mechanism(conditionsTup[1])
        return results
    except:
        results='failed'
        return results

if len(conditionsTups)>8:
    cores=8
else:
    cores=len(conditionsTups)
    
pool = ThreadPool(cores) 
results = pool.map(solver,conditionsTups)

