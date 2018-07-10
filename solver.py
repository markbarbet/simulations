# -*- coding: utf-8 -*-
"""
Created on Tue Jun 05 14:59:37 2018

@author: HP USER
"""

import cantera as ct
import simulations as sim

def JSR(filenames,Temps,pressures,resTime,conditions,sens,volume,observables,fuel):
    expConditions=[]

    
    for T in Temps:
        for filename in filenames:
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
        #results[-1].add_fuel(fuel)
        #results[-1].assign_phi(phi)
        results[-1].add_residence_time(i[3])
        #print(results[-1].final_pressure)
    return results