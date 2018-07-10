# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:02:35 2017

@author: Carly LaGrotta
"""
import itertools
import numpy as np
import cantera as ct
import pandas as pd
import simulations
        

def ShockTube(ctiFile ,pressure,temperature,conditions,initialTime,finalTime,thermalBoundary,observables=[],physical_params=[], kinetic_sens=0,physical_sens=0):
    #gas = ct.Solution('AramcoMech2.0.cti')
    # 'GRI30-1999.cti'
    gas = ct.Solution(ctiFile)
    gas.TPX = temperature, pressure*101325, conditions
    physical_params = ['T','P','X']
    
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
                
                if(newcounter%gas.n_reactions == 0):
                    dfs[observables.index(observable)] = dfs[observables.index(observable)].append(((pd.DataFrame(tempArray[observables.index(observable)])).transpose()),ignore_index=True)
                
        counter +=1
    if kinetic_sens == 1 and bool(observables):
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
                    

                
                newTime = 0 
                newCounter = 0
                while newTime < finalTime:
                    newTime = sim2.step()
                    state = np.hstack([newTime, shockTube2.thermo.P, shockTube2.mass, 
                    shockTube2.volume, shockTube2.T, shockTube2.thermo.X])
                    timeHistorytest.loc[newCounter] = state
                    newCounter +=1

                newTimeArray = timeHistorytest['time'] 
                tempForInterp = timeHistorytest[observables]
            
                tempForInterplstT = [tempForInterp.ix[:,x].values for x in xrange(tempForInterp.shape[1])]
                
                interpolatedData = [np.interp(timeHistory['time'].values,newTimeArray.values,tempForInterplstT[x]) for x in xrange(len(tempForInterplstT))]
                
                interpolatedData = [pd.DataFrame(interpolatedData[x]) for x in xrange(len(interpolatedData))]
                interpolatedData = pd.concat(interpolatedData, axis=1,ignore_index=True)
                
                

                tempT = interpolatedData.applymap(np.log10)
                tempT.columns = observables
                tempT = (originalPsens.subtract(tempT))/np.log10(dk)
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
                    
                
                newTime2 = 0
                newCounter2 = 0
                while newTime2 < finalTime:
                    newTime2 = sim3.step()
                    state = np.hstack([newTime2, shockTube3.thermo.P, shockTube3.mass, 
                    shockTube3.volume, shockTube3.T, shockTube3.thermo.X])
                    timeHistorytest2.loc[newCounter2] = state
                    newCounter2 +=1
                    
                newTimeArray2 = timeHistorytest2['time']

                tempForInterp = timeHistorytest2[observables]
                tempForInterplstP = [tempForInterp.ix[:,x].values for x in xrange(tempForInterp.shape[1])]
                                     
                interpolatedData2 = [np.interp(timeHistory['time'].values,newTimeArray2.values,tempForInterplstP[x]) for x in xrange(len(tempForInterplstP))]
                                     
                interpolatedData2 = [pd.DataFrame(interpolatedData2[x]) for x in xrange(len(interpolatedData2))]
                interpolatedData2 = pd.concat(interpolatedData2,axis=1,ignore_index=True)
                
                tempP = interpolatedData2.applymap(np.log10)
                tempP.columns = observables
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
                

                newTime3 = 0
                newCounter3 = 0
                while newTime3 < finalTime:
                    newTime3 = sim4.step()                    
                    state = np.hstack([newTime3, shockTube4.thermo.P, shockTube4.mass, 
                    shockTube4.volume, shockTube4.T, shockTube4.thermo.X])
                    
                    timeHistorytest3.loc[newCounter3] = state    
                    newCounter3 +=1
                newTimeArray3 = timeHistorytest3['time']
                tempForInterp = timeHistorytest3[observables]
                tempForInterplstX = [tempForInterp.ix[:,x].values for x in xrange(tempForInterp.shape[1])]
                
                interpolatedData3 = [np.interp(timeHistory['time'].values,newTimeArray3.values,tempForInterplstX[x]) for x in xrange(len(tempForInterplstX))]
                interpolatedData3 = [pd.DataFrame(interpolatedData3[x]) for x in xrange(len(interpolatedData3))]
                interpolatedData3 = pd.concat(interpolatedData3,axis=1,ignore_index = True)
                tempX = interpolatedData3.applymap(np.log10)
                tempX.columns = observables
                tempX = (originalPsens.subtract(tempX))/np.log10(dk)
                tempXlst = [tempX.ix[:,idx] for idx in xrange(tempX.shape[1])]

                                

    if physical_sens == 1 and bool(observables):         

        if 'T' in physical_params and 'P' in physical_params and 'X' in physical_params:
            t = [tempTlst,tempPlst,tempXlst]
            psensIndex = [timeHistory['time'].as_matrix(),['T','P','X'],observables]
            psensdfs = [pd.concat([t[0][x],t[1][x],t[2][x]],ignore_index = True , axis = 1) for x in xrange(len(tempXlst))]
            numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
            pS = np.dstack(numpyMatrixspsens)
            
        elif 'T' in physical_params and 'P' in physical_params:
            t = [tempTlst,tempPlst]
            psensIndex = [timeHistory['time'].as_matrix(),['T','P'],observables]
            psensdfs = [pd.concat([t[0][x],t[1][x]],ignore_index=True,axis = 1) for x in xrange(len(tempTlst))]                
            numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
            pS = np.dstack(numpyMatrixspsens)
            
        elif 'T'in physical_params and 'X' in physical_params:
            t = [tempTlst,tempXlst]
            psensIndex = [timeHistory['time'].as_matrix(),['T','X'],observables]
            psensdfs = [pd.concat([t[0][x],t[1][x]], ignore_index = True, axis = 1) for x in xrange(len(tempTlst))]
            numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
            pS = np.dstack(numpyMatrixspsens)
            
        elif 'P'in physical_params and 'X' in physical_params:
            t = [tempPlst,tempXlst]
            psensIndex = [timeHistory['time'].as_matrix(),['P','X'],observables]
            psensdfs = [pd.concat([t[0][x],t[1][x]], ignore_index = True, axis = 1) for x in xrange(len(tempPlst))]
            numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
            pS = np.dstack(numpyMatrixspsens)
            
        elif 'T' in physical_params:
            t = [tempTlst]
            psensIndex = [timeHistory['time'].as_matrix(),['T'],observables]
            numpyMatrixspsens = [tempTlst[dataframe].as_matrix() for dataframe in xrange(len(tempTlst))]
            pS = np.dstack(numpyMatrixspsens)
            
        elif 'P' in physical_params:
            t = [tempPlst]
            psensIndex = [timeHistory['time'].as_matrix(),['P'],observables]
            numpyMatrixspsens = [tempPlst[dataframe].as_matrix() for dataframe in xrange(len(tempPlst))]
            pS = np.dstack(numpyMatrixspsens)
            
        elif 'X' in physical_params:
            t = [tempXlst]
            psensIndex = [timeHistory['time'].as_matrix(),['X'],observables]
            numpyMatrixspsens = [tempXlst[dataframe].as_matrix() for dataframe in xrange(len(tempXlst))]
            pS = np.dstack(numpyMatrixspsens)

        

                    
    #print(np.shape(pS))   
    if kinetic_sens==1 and bool(observables) and physical_sens==0:  
        results = simulations.model_data('Shock-Tube',kinetic_sens = S, Solution = timeHistory, Index = ksensIndex)
        return results 
        
    if kinetic_sens==0 and bool(physical_params) and physical_sens==1:
        results = simulations.model_data('Shock-Tube',physical_sens = pS, Solution = timeHistory, pIndex = psensIndex)
        return results 

    if physical_sens==1 and bool(observables)and kinetic_sens==1:
        results = simulations.model_data('Shock-Tube',kinetic_sens = S, physical_sens = pS, Solution = timeHistory, Index = ksensIndex, pIndex = psensIndex)
        return results 
      
    if physical_sens == 0 and kinetic_sens == 0:
        results = simulations.model_data('Shock-Tube',Solution = timeHistory)
        return results