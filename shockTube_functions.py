# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:02:35 2017

@author: Carly LaGrotta
"""
import sys, itertools
import numpy as np
import cantera as ct
import pandas as pd
import yaml


        
class shock_tube_model_data:
    def __init__(self,kinetic_sens=[],physical_sens=[],Solution=[],Index=[],pIndex=[]):
        self.k_sens = kinetic_sens
        self.p_sens = physical_sens
        self.solution = Solution
        self.Index = Index 
        self.pIndex = pIndex 
    
def ShockTube(ctiFile ,pressure,temperature,conditions,initialTime,finalTime,thermalBoundary,observables=[],physical_params=[], kinetic_sens=0,physical_sens=0):
    #gas = ct.Solution('AramcoMech2.0.cti')
    # 'GRI30-1999.cti'
    gas = ct.Solution(ctiFile)
    gas.TPX = temperature, pressure*101325, conditions
    
    if thermalBoundary == 'adiabatic': 
        shockTube = ct.IdealGasReactor(gas,name = 'R1',energy= 'on')
    elif thermalBoundary == 'isothermal': 
        shockTube = ct.IdealGasReactor(gas,name = 'R1', energy= 'off')
    else:
        raise Exception('Please enter adiabatic or isothermal for the thermal boundary layer')
    
    
    sim=ct.ReactorNet([shockTube])
    

 
    if kinetic_sens==1 and bool(observables)==False:
        raise Exception('Please supply a non-empty list of observables for sensitivity analysis or set kinetic_sens=0')
        
        
    if physical_sens==1 and bool(observables)==False:
        raise Exception('Please supply a non-empty list of observables for sensitivity analysis or set physical_sens=0')
    
    
    if kinetic_sens==1 and bool(observables):
        [shockTube.add_sensitivity_reaction(i) for i in xrange(gas.n_reactions)]
        dfs = [pd.DataFrame() for x in xrange(len(observables))]
        tempArray = [np.zeros(gas.n_reactions) for x in xrange(len(observables))]
                
    
    if physical_sens==1 and bool(observables):
        baseConditions = gas.TPX
                      
        
    
    columnNames = [shockTube.component_name(item) for item in range(shockTube.n_vars)]               
    columnNames = ['time']+['pressure'] + columnNames
    timeHistory = pd.DataFrame(columns = columnNames)
    timeHistorytest= pd.DataFrame(columns = columnNames)
    timeHistorytest2= pd.DataFrame(columns = columnNames)
    timeHistorytest3 = pd.DataFrame(columns = columnNames)
    
   
    

    
    t=initialTime
    counter = 0
    while t < finalTime:
        
        t = sim.step() 
        
        state = np.hstack([t, shockTube.thermo.P, shockTube.mass, 
                   shockTube.volume, shockTube.T, shockTube.thermo.X])
        timeHistory.loc[counter] = state

        if kinetic_sens == 1 and bool(observables):

            newcounter = 0
            for observable,reaction in itertools.product(observables,xrange(gas.n_reactions)):

                tempArray[observables.index(observable)][reaction] = sim.sensitivity(observable, reaction)
               
                newcounter +=1
                
                #print('this is new counter' , newcounter)
                
                if(newcounter%gas.n_reactions == 0):
             
                    print(observable , newcounter)
                    dfs[observables.index(observable)] = dfs[observables.index(observable)].append(((pd.DataFrame(tempArray[observables.index(observable)])).transpose()),ignore_index=True)
                
        counter +=1

    ksensIndex = [timeHistory['time'].as_matrix(),gas.reaction_equations(),observables]
   
    
    numpyMatrixsksens = [dfs[dataframe].as_matrix() for dataframe in xrange(len(dfs))]
    S = np.dstack(numpyMatrixsksens)
    print(np.shape(S))
        
    if physical_sens == 1 and bool(observables):
        dk=.01
        originalPsens = (timeHistory[observables]).applymap(np.log10)
        
        
        
        for numOfPhyParams in xrange(len(physical_params)):
            if physical_params[numOfPhyParams] == 'T':
                gas2 = ct.Solution(ctiFile)
                gas2.TPX = baseConditions[0]+dk,baseConditions[1],baseConditions[2]
                if thermalBoundary == 'adiabatic': 
                    shockTube2 = ct.IdealGasReactor(gas2,name = 'R1',energy= 'on')
                    sim2=ct.ReactorNet([shockTube2])
                if thermalBoundary == 'isothermal':
                    shockTube2 = ct.IdealGasReactor(gas2,name = 'R1',energy= 'off')
                    sim2=ct.ReactorNet([shockTube2])
                    

                for idx, newTime in enumerate(timeHistory['time'].as_matrix()):
                    sim2.advance(newTime)
                    state = np.hstack([newTime, shockTube2.thermo.P, shockTube2.mass, 
                    shockTube2.volume, shockTube2.T, shockTube2.thermo.X])
                    timeHistorytest.loc[idx] = state
                tempT = (timeHistorytest[observables]).applymap(np.log10)
                tempT = (originalPsens.subtract(tempT))/np.log10(dk)
                print(tempT)
                tempTlst = [tempT.ix[:,idx] for idx in xrange(tempT.shape[1])]
                
                 
                
            if physical_params[numOfPhyParams] == 'P':
                gas3 = ct.Solution(ctiFile)
                gas3.TPX = baseConditions[0],baseConditions[1]+dk,baseConditions[2]
                if thermalBoundary == 'adiabatic':
                    shockTube3 = ct.IdealGasReactor(gas3,name = 'R1',energy = 'on')
                    sim3 = ct.ReactorNet([shockTube3])
                if thermalBoundary =='isothermal':
                    shockTube3 = ct.IdealGasReactor(gas3,name = 'R1', energy = 'off')
                    sim3 = ct.ReactorNet([shockTube3])
                    
                for idx2, newTime2 in enumerate(timeHistory['time'].as_matrix()):
                    sim3.advance(newTime)
                    state = np.hstack([newTime2, shockTube3.thermo.P, shockTube3.mass, 
                    shockTube3.volume, shockTube3.T, shockTube3.thermo.X])
                    timeHistorytest2.loc[idx2] = state

                tempP = (timeHistorytest2[observables]).applymap(np.log10)
                tempP = (originalPsens.subtract(tempP))/np.log10(dk)
                tempPlst = [tempP.ix[:,idx] for idx in xrange(tempP.shape[1])]
                
                
                            
            if physical_params[numOfPhyParams] == 'X':
                gas4 = ct.Solution(ctiFile)
                gas4.TPX = baseConditions[0],baseConditions[1],baseConditions[2]+dk
                if thermalBoundary =='adiabatic':
                    shockTube4 = ct.IdealGasReactor(gas4,name = 'R1',energy= 'on')
                    sim4 = ct.ReactorNet([shockTube4])
                if thermalBoundary == 'isothermal':
                    shockTube4 = ct.IdealGasReactor(gas4 ,name = 'R1', energy = 'off')
                    sim4 = ct.ReactorNet([shockTube4]) 
                

                for idx3, newTime3 in enumerate(timeHistory['time'].as_matrix()):
                    sim4.advance(newTime3)
                    state = np.hstack([newTime3, shockTube4.thermo.P, shockTube4.mass, 
                    shockTube4.volume, shockTube4.T, shockTube4.thermo.X])
                    timeHistorytest3.loc[idx3] = state      
      
                tempX = (timeHistorytest3[observables]).applymap(np.log10)
                tempX = (originalPsens.subtract(tempX))/np.log10(dk)
                tempXlst = [tempX.ix[:,idx] for idx in xrange(tempX.shape[1])]                    

    
 
        

    if "T" and "P" and "X" in physical_params:
        t = [tempTlst,tempPlst,tempXlst]
        psensIndex = [timeHistory['time'].as_matrix(),['T','P','X'],observables]
        psensdfs = [pd.concat([t[0][x],t[1][x],t[2][x]],ignore_index = True , axis = 1) for x in xrange(len(tempXlst))]
        numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
        pS = np.dstack(numpyMatrixspsens)
        
    elif "T" and "P" in physical_params:
        t = [tempTlst,tempPlst]
        psensIndex = [timeHistory['time'].as_matrix(),['T','P'],observables]
        psensdfs = [pd.concat([t[0][x],t[1][x]],ignore_index=True,axis = 1) for x in xrange(len(tempTlst))]                
        numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
        pS = np.dstack(numpyMatrixspsens)
        print(np.shape(pS))
        
        
    elif "T" and "X" in physical_params:
        t = [tempTlst,tempXlst]
        psensIndex = [timeHistory['time'].as_matrix(),['T','X'],observables]
        psensdfs = [pd.concat([t[0][x],t[1][x]], ignore_index = True, axis = 1) for x in xrange(len(tempTlst))]
        numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
        pS = np.dstack(numpyMatrixspsens)
        
    elif "P" and "X" in physical_params:
        t = [tempPlst,tempXlst]
        psensIndex = [timeHistory['time'].as_matrix(),['P','X'],observables]
        psensdfs = [pd.concat([t[0][x],t[1][x]], ignore_index = True, axis = 1) for x in xrange(len(tempPlst))]
        numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
        pS = np.dstack(numpyMatrixspsens)
        
    elif "T" in physical_params:
        t = [tempTlst]
        psensIndex = [timeHistory['time'].as_matrix(),['T'],observables]
        numpyMatrixspsens = [tempTlst[dataframe].as_matrix() for dataframe in xrange(len(tempTlst))]
        pS = np.dstack(numpyMatrixspsens)
        
    elif "P" in physical_params:
        t = [tempPlst]
        psensIndex = [timeHistory['time'].as_matrix(),['P'],observables]
        numpyMatrixspsens = [tempPlst[dataframe].as_matrix() for dataframe in xrange(len(tempPlst))]
        pS = np.dstack(numpyMatrixspsens)
        
    elif "X" in physical_params:
        t = [tempXlst]
        psensIndex = [timeHistory['time'].as_matrix(),['X'],observables]
        numpyMatrixspsens = [tempXlst[dataframe].as_matrix() for dataframe in xrange(len(tempPlst))]
        pS = np.dstack(numpyMatrixspsens)

        

                
 
    
    
                 
    
        
 
        
    if kinetic_sens==1 and bool(observables) and physical_sens==0:  
        results = shock_tube_model_data(kinetic_sens = S, Solution = timeHistory, Index = ksensIndex)
        return results 
        
    if kinetic_sens==0 and bool(physical_params) and physical_sens==1:
        results = shock_tube_model_data(physical_sens = pS, solution = timeHistory, pIndex = psensIndex)
        return results 

    if physical_sens==1 and bool(observables)and kinetic_sens==1:
        results = shock_tube_model_data(kinetic_sens = S, physical_sens = pS, Solution = timeHistory, Index = ksensIndex, pIndex = psensIndex)
        return results 
      
    if physical_sens == 0 and kinetic_sens == 0:
        results = shock_tube_model_data(Solution = timeHistory)
        return results