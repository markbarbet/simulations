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
        

def ShockTube(ctiFile,speciesNames,pressure,temperature,conditions,initialTime,finalTime,thermalBoundary,observables=[],physical_params=[], kinetic_sens=0,physical_sens=0,reactorType='cv'):
    #gas = ct.Solution('AramcoMech2.0.cti')
    # 'GRI30-1999.cti'
    if 'AR' in speciesNames:
        addBackin = 'AR'
        speciesNames.remove('AR')
    if 'Ar' in speciesNames:
        addBackin = 'Ar'
        speciesNames.remove('Ar')
    if 'HE' in speciesNames:
        addBackin = 'HE'
        speciesNames.remove('HE')
    if 'He' in speciesNames:
        addBackin = 'He'
        speciesNames.remove('He')
            
    
    
    gas = ct.Solution(ctiFile)
    gas.TPX = temperature, pressure*101325, conditions
    physicalParamsSpecies = ['X'+ species for species in speciesNames]
    
    physical_params = ['P','T','X']
    
    if thermalBoundary == 'adiabatic':
        if reactorType=='cv':
            shockTube = ct.IdealGasReactor(gas,name = 'R1',energy= 'on')
        elif reactorType=='cp':
            shockTube = ct.IdealGasConstPressureReactor(gas,name = 'R1', energy='on')
    elif thermalBoundary == 'isothermal': 
        if reactorType=='cv':
            shockTube = ct.IdealGasReactor(gas,name = 'R1', energy= 'off')
        elif reactorType=='cp':
            shockTube = ct.IdealGasConstPressureReactor(gas,name = 'R1', energy = 'off')
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
        originalConditions = conditions 
                      
        
    
    columnNames = [shockTube.component_name(item) for item in xrange(shockTube.n_vars)]               
    columnNames = ['time']+['pressure'] + columnNames
    timeHistory = pd.DataFrame(columns = columnNames)
    timeHistorytest= pd.DataFrame(columns = columnNames)
    timeHistorytest2= pd.DataFrame(columns = columnNames)
    timeHistorytest4 = [pd.DataFrame(columns = columnNames) for item in xrange(len(physicalParamsSpecies))]
    
   
    t=initialTime
    counter = 0
    #commenting out volume
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
        numpyMatrixsBsens = [numpyMatrixsksens[x]*(np.log(timeHistory['temperature'].as_matrix().flatten()))[:,np.newaxis] for x in xrange(len(numpyMatrixsksens))]
        numpyMatrixsEsens = [numpyMatrixsksens[x]*(np.true_divide(-1,timeHistory['temperature'].as_matrix().flatten()))[:,np.newaxis] for x in xrange(len(numpyMatrixsksens))]
        S = np.dstack(numpyMatrixsksens)
        Sa = S
        Sb = np.dstack(numpyMatrixsBsens)
        Se = np.dstack(numpyMatrixsEsens)                             
        S = np.dstack(numpyMatrixsksens)
       
       
    if physical_sens == 1 and bool(observables):
        dk=.01
        originalPsens = (timeHistory[observables]).applymap(np.log)
                        
        for numOfPhyParams in xrange(len(physical_params)):
            if physical_params[numOfPhyParams] == 'T':
                gas2 = ct.Solution(ctiFile)
                gas2.TPX = baseConditions[0]+dk,baseConditions[1],baseConditions[2]
                if thermalBoundary == 'adiabatic': 
                    if reactorType=='cv':
                        shockTube2 = ct.IdealGasReactor(gas2,name = 'R1',energy= 'on')
                    elif reactorType=='cp':
                        shockTube2 = ct.IdealGasConstPressureReactor(gas2,name = 'R1', energy='on')
                    sim2=ct.ReactorNet([shockTube2])
                if thermalBoundary == 'isothermal':
                    if reactorType=='cv':
                        shockTube2 = ct.IdealGasReactor(gas2,name = 'R1',energy= 'off')
                    elif reactorType=='cp':
                        shockTube2 = ct.IdealGasConstPressureReactor(gas2,name = 'R1', energy='off')
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
                temperatureForMappingPhysicalSensT = timeHistorytest['temperature'].values
                pressureForMappingPhysicalSensT = timeHistorytest['pressure'].values


                tempForInterp = timeHistorytest[observables]
            
                tempForInterplstT = [tempForInterp.ix[:,x].values for x in xrange(tempForInterp.shape[1])]
                
                interpolatedData = [np.interp(timeHistory['time'].values,newTimeArray.values,tempForInterplstT[x]) for x in xrange(len(tempForInterplstT))]
                interpolatedTemperatureForMappingPhysicalSensT = [np.interp(timeHistory['time'].values,newTimeArray.values,temperatureForMappingPhysicalSensT)]
                interpolatedPressureForMappingPhysicalSensT = [np.interp(timeHistory['time'].values,newTimeArray.values,pressureForMappingPhysicalSensT)]
                interpolatedData = [pd.DataFrame(interpolatedData[x]) for x in xrange(len(interpolatedData))]
                interpolatedData = pd.concat(interpolatedData, axis=1,ignore_index=True)
                
                
                concentrationOfAbsorbanceObservablesForSensT = interpolatedData
                concentrationOfAbsorbanceObservablesForSensT.columns = observables
                concentrationOfAbsorbanceObservablesForSensT = [concentrationOfAbsorbanceObservablesForSensT]
                
                

                tempT = interpolatedData.applymap(np.log)
                tempT.columns = observables
                tempT = (originalPsens.subtract(tempT))/np.log(dk)
                tempTlst = [tempT.ix[:,idx] for idx in xrange(tempT.shape[1])]
                
                
                
                
                 
                
            if physical_params[numOfPhyParams] == 'P':
                gas3 = ct.Solution(ctiFile)
                gas3.TPX = baseConditions[0],baseConditions[1]+dk,baseConditions[2]
                if thermalBoundary == 'adiabatic':
                    if reactorType=='cv':
                        shockTube3 = ct.IdealGasReactor(gas3,name = 'R1',energy = 'on')
                    elif reactorType=='cp':
                        shockTube3 = ct.IdealGasConstPressureReactor(gas3,name = 'R1', energy='on')
                    sim3 = ct.ReactorNet([shockTube3])
                if thermalBoundary =='isothermal':
                    if reactorType=='cv':                        
                        shockTube3 = ct.IdealGasReactor(gas3,name = 'R1', energy = 'off')
                    elif reactorType=='cp':
                        shockTube3 = ct.IdealGasConstPressureReactor(gas3,name = 'R1', energy='off')
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
                temperatureForMappingPhysicalSensP = timeHistorytest2['temperature'].values
                pressureForMappingPhysicalSensP = timeHistorytest2['pressure'].values

                tempForInterp = timeHistorytest2[observables]
                tempForInterplstP = [tempForInterp.ix[:,x].values for x in xrange(tempForInterp.shape[1])]
                                     
                interpolatedData2 = [np.interp(timeHistory['time'].values,newTimeArray2.values,tempForInterplstP[x]) for x in xrange(len(tempForInterplstP))]
                interpolatedTemperatureForMappingPhysicalSensP = [np.interp(timeHistory['time'].values,newTimeArray2.values,temperatureForMappingPhysicalSensP)]
                
                
                interpolatedPressureForMappingPhysicalSensP = [np.interp(timeHistory['time'].values,newTimeArray2.values,pressureForMappingPhysicalSensP)]                     
                interpolatedData2 = [pd.DataFrame(interpolatedData2[x]) for x in xrange(len(interpolatedData2))]
                interpolatedData2 = pd.concat(interpolatedData2,axis=1,ignore_index=True)
                
                concentrationOfAbsorbanceObservablesForSensP = interpolatedData2
                concentrationOfAbsorbanceObservablesForSensP.columns = observables
                concentrationOfAbsorbanceObservablesForSensP = [concentrationOfAbsorbanceObservablesForSensP]
                
                
                tempP = interpolatedData2.applymap(np.log)
                tempP.columns = observables
                tempP = (originalPsens.subtract(tempP))/np.log(dk)
                tempPlst = [tempP.ix[:,idx] for idx in xrange(tempP.shape[1])]
                                     
        
                
                            
            if physical_params[numOfPhyParams] == 'X':
                for simulationNumber ,  speciesName  in enumerate(speciesNames):
                    
                    
                    newConditions = originalConditions
                    gas4 = ct.Solution(ctiFile)
                    originalValue = conditions[speciesName]
                    newValue = originalValue + dk
                    newConditions.update({speciesName:newValue})
                    gas4.TPX = baseConditions[0],baseConditions[1], newConditions
                   
                    if thermalBoundary =='adiabatic':
                        if reactorType=='cv':
                            shockTube4 = ct.IdealGasReactor(gas4,name = 'R1',energy= 'on')
                        elif reactorType=='cp':
                            shockTube4 = ct.IdealGasConstPressureReactor(gas4,name = 'R1', energy='on')
                        sim4 = ct.ReactorNet([shockTube4])
                    if thermalBoundary == 'isothermal':
                        if reactorType=='cv':
                            shockTube4 = ct.IdealGasReactor(gas4 ,name = 'R1', energy = 'off')
                        elif reactorType=='cp':
                            shockTube4 = ct.IdealGasConstPressureReactor(gas4,name = 'R1', energy='off')
                        sim4 = ct.ReactorNet([shockTube4]) 
                    
    
                    newTime3 = 0
                    newCounter3 = 0
                    while newTime3 < finalTime:
                        newTime3 = sim4.step()                    
                        state = np.hstack([newTime3, shockTube4.thermo.P, shockTube4.mass, 
                        shockTube4.volume, shockTube4.T, shockTube4.thermo.X])
                        
                        timeHistorytest4[simulationNumber].loc[newCounter3] = state    
                        newCounter3 +=1
                    
                    newConditions.update({speciesName:originalValue})
                
                
                    
                newTimeArrays3 = [timeHistorytest4[simulationNumber]['time'] for simulationNumber in xrange(len(timeHistorytest4))]
                temperatureForMappingPhysicalSensX = [timeHistorytest4[simulationNumber]['temperature'] for simulationNumber in xrange(len(timeHistorytest4))]
                pressureForMappingPhysicalSensX = [timeHistorytest4[simulationNumber]['pressure'] for simulationNumber in xrange(len(timeHistorytest4))]                
                tempForInterps = [timeHistorytest4[simulationNumber][observables] for simulationNumber in xrange(len(timeHistorytest4))]
                
                
                tempForInterplstXs = [[] for x in xrange(len(timeHistorytest4))]
                for simulationNumber in xrange(len(timeHistorytest4)):
                    for  x in xrange(tempForInterps[simulationNumber].shape[1]):
                        tempForInterplstXs[simulationNumber].append(tempForInterps[simulationNumber].ix[:,x].values)
                        

                
                interpolatedData3 = [[] for x in xrange(len(timeHistorytest4))]
                interpolatedTemperatureForMappingPhysicalSensX = []
                interpolatedPressureForMappingPhysicalSensX = [] 
               
                                     
                for simulationNumber in xrange(len(timeHistorytest4)):
                    for x in xrange(tempForInterps[simulationNumber].shape[1]):
                        interpolatedData3[simulationNumber].append(np.interp(timeHistory['time'].values,newTimeArrays3[simulationNumber].values,tempForInterplstXs[simulationNumber][x]))
                for simulationNumber in xrange(len(timeHistorytest4)):
                    interpolatedTemperatureForMappingPhysicalSensX.append(np.interp(timeHistory['time'].values,newTimeArrays3[simulationNumber].values,temperatureForMappingPhysicalSensX[simulationNumber]))
                
                for simulationNumber in xrange(len(timeHistorytest4)):
                    interpolatedPressureForMappingPhysicalSensX.append(np.interp(timeHistory['time'].values,newTimeArrays3[simulationNumber].values,pressureForMappingPhysicalSensX[simulationNumber]))                        
                        
                interpolatedDataFrames = [[] for x in xrange(len(timeHistorytest4))]

                #print(interpolatedData3)
                for simulationNumber in xrange(len(timeHistorytest4)):
                    for x in xrange(tempForInterps[simulationNumber].shape[1]):
                        interpolatedDataFrames[simulationNumber].append(pd.DataFrame(interpolatedData3[simulationNumber][x]))
                
                interpolatedDataFrames = [pd.concat(interpolatedDataFrames[simulationNumber],axis=1,ignore_index=True) for simulationNumber in xrange(len(timeHistorytest4))]
                concentrationOfAbsorbanceObservablesForSensX = interpolatedDataFrames
                
                for simulationNumber in xrange(len(timeHistorytest4)):
                    concentrationOfAbsorbanceObservablesForSensX[simulationNumber].columns = observables
                    
                    
                interpolatedDataFrames = [interpolatedDataFrames[simulationNumber].applymap(np.log) for simulationNumber in xrange(len(timeHistorytest4))]
                for simulationNumber in xrange(len(timeHistorytest4)):
                    interpolatedDataFrames[simulationNumber].columns = observables
                interpolatedDataFrames = [(originalPsens.subtract(interpolatedDataFrames[simulationNumber]))/np.log(dk) for simulationNumber in xrange(len(timeHistorytest4))]
                                          
                tempXlst = [[] for x in xrange(len(timeHistorytest4))]
                for simulationNumber in xrange(len(timeHistorytest4)):
                    for x in xrange(interpolatedDataFrames[simulationNumber].shape[1]):
                        tempXlst[simulationNumber].append(interpolatedDataFrames[simulationNumber].ix[:,x])
                        
              
                tempXlsts = [[] for x in xrange(len(observables))]
                             
                for lengthOfList in xrange(len(observables)):
                    for dfInList in xrange(len(timeHistorytest4)):
                        tempXlsts[lengthOfList].append(tempXlst[dfInList][lengthOfList])
                        
                tempXlsts = [pd.concat(tempXlsts[simulationNumber],axis=1,ignore_index=True) for simulationNumber in xrange(len(tempXlsts))]
                
 



                                

    if physical_sens == 1 and bool(observables):         
        speciesNames.append(addBackin)
        if 'T' in physical_params and 'P' in physical_params and 'X' in physical_params:
            t = [tempTlst,tempPlst,tempXlsts]
            psensIndex = [timeHistory['time'].as_matrix(),['T','P'] + physicalParamsSpecies ,observables]
            psensdfs = [pd.concat([t[0][x],t[1][x],t[2][x]],ignore_index = True , axis = 1) for x in xrange(len(tempXlsts))]
            numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
            interpolatedTemperatureForMappingPhysicalSens = [interpolatedTemperatureForMappingPhysicalSensT + interpolatedTemperatureForMappingPhysicalSensP + interpolatedTemperatureForMappingPhysicalSensX]
            interpolatedPressureForMappingPhysicalSens = [interpolatedPressureForMappingPhysicalSensT + interpolatedPressureForMappingPhysicalSensP + interpolatedPressureForMappingPhysicalSensX] 
            concentrationOfAbsorbanceObservablesForSens = [concentrationOfAbsorbanceObservablesForSensT+concentrationOfAbsorbanceObservablesForSensP+concentrationOfAbsorbanceObservablesForSensX]

            pS = np.dstack(numpyMatrixspsens)
#            
        elif 'T' in physical_params and 'P' in physical_params:
            t = [tempTlst,tempPlst]
            psensIndex = [timeHistory['time'].as_matrix(),['T','P'],observables]
            psensdfs = [pd.concat([t[0][x],t[1][x]],ignore_index=True,axis = 1) for x in xrange(len(tempTlst))]                
            numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
            interpolatedTemperatureForMappingPhysicalSens = [interpolatedTemperatureForMappingPhysicalSensT + interpolatedTemperatureForMappingPhysicalSensP]
            interpolatedPressureForMappingPhysicalSens = [interpolatedPressureForMappingPhysicalSensT + interpolatedPressureForMappingPhysicalSensP ]     
            concentrationOfAbsorbanceObservablesForSens = [concentrationOfAbsorbanceObservablesForSensT+concentrationOfAbsorbanceObservablesForSensP]

            pS = np.dstack(numpyMatrixspsens)
            
        elif 'T'in physical_params and 'X' in physical_params:
            t = [tempTlst,tempXlsts]
            psensIndex = [timeHistory['time'].as_matrix(),['T'] + physicalParamsSpecies ,observables]
            psensdfs = [pd.concat([t[0][x],t[1][x]], ignore_index = True, axis = 1) for x in xrange(len(tempTlst))]
            numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
            interpolatedTemperatureForMappingPhysicalSens = [interpolatedTemperatureForMappingPhysicalSensT + interpolatedTemperatureForMappingPhysicalSensX]
            interpolatedPressureForMappingPhysicalSens = [interpolatedPressureForMappingPhysicalSensT +interpolatedPressureForMappingPhysicalSensX] 
            concentrationOfAbsorbanceObservablesForSens = [concentrationOfAbsorbanceObservablesForSensT+concentrationOfAbsorbanceObservablesForSensX]

            pS = np.dstack(numpyMatrixspsens)
#            
        elif 'P'in physical_params and 'X' in physical_params:
            t = [tempPlst,tempXlsts]
            psensIndex = [timeHistory['time'].as_matrix(),['P'] + physicalParamsSpecies ,observables]
            psensdfs = [pd.concat([t[0][x],t[1][x]], ignore_index = True, axis = 1) for x in xrange(len(tempPlst))]
            numpyMatrixspsens = [psensdfs[dataframe].as_matrix() for dataframe in xrange(len(psensdfs))]
            interpolatedTemperatureForMappingPhysicalSens = [interpolatedTemperatureForMappingPhysicalSensP + interpolatedTemperatureForMappingPhysicalSensX]
            interpolatedPressureForMappingPhysicalSens = [interpolatedPressureForMappingPhysicalSensP + interpolatedPressureForMappingPhysicalSensX]
            concentrationOfAbsorbanceObservablesForSens = [concentrationOfAbsorbanceObservablesForSensP+concentrationOfAbsorbanceObservablesForSensX]

            pS = np.dstack(numpyMatrixspsens)
            
        elif 'T' in physical_params:
            t = [tempTlst]
            psensIndex = [timeHistory['time'].as_matrix(),['T'],observables]
            numpyMatrixspsens = [tempTlst[dataframe].as_matrix() for dataframe in xrange(len(tempTlst))]
            interpolatedTemperatureForMappingPhysicalSens = [interpolatedTemperatureForMappingPhysicalSensT]
            interpolatedPressureForMappingPhysicalSens = [interpolatedPressureForMappingPhysicalSensT ]
            concentrationOfAbsorbanceObservablesForSens = [concentrationOfAbsorbanceObservablesForSensT]

            pS = np.dstack(numpyMatrixspsens)
            
        elif 'P' in physical_params:
            t = [tempPlst]
            psensIndex = [timeHistory['time'].as_matrix(),['P'],observables]
            numpyMatrixspsens = [tempPlst[dataframe].as_matrix() for dataframe in xrange(len(tempPlst))]
            interpolatedTemperatureForMappingPhysicalSens = [interpolatedTemperatureForMappingPhysicalSensP ]
            interpolatedPressureForMappingPhysicalSens = [interpolatedPressureForMappingPhysicalSensP ]
            concentrationOfAbsorbanceObservablesForSens = [concentrationOfAbsorbanceObservablesForSensP]

            pS = np.dstack(numpyMatrixspsens)
            
        elif 'X' in physical_params:
            t = [tempXlsts]
            psensIndex = [timeHistory['time'].as_matrix(),[physicalParamsSpecies],observables]
            numpyMatrixspsens = [tempXlsts[dataframe].as_matrix() for dataframe in xrange(len(tempXlsts))]
            interpolatedTemperatureForMappingPhysicalSens = [interpolatedTemperatureForMappingPhysicalSensX]
            interpolatedPressureForMappingPhysicalSens = [interpolatedPressureForMappingPhysicalSensX]
            concentrationOfAbsorbanceObservablesForSens = [concentrationOfAbsorbanceObservablesForSensX]

            pS = np.dstack(numpyMatrixspsens)
                   
 
    if kinetic_sens==1 and bool(observables) and physical_sens==0:  
        results = simulations.model_data('Shock-Tube',kinetic_sens = S, Solution = timeHistory, Index = ksensIndex)
        results.assign_ksens_mappings(Sa ,Sb, Se)
        return results 
        
    if kinetic_sens==0 and bool(physical_params) and physical_sens==1:
        results = simulations.model_data('Shock-Tube',physical_sens = pS, Solution = timeHistory, pIndex = psensIndex)
        results.add_interpolated_temps_for_physical_sens(interpolatedTemperatureForMappingPhysicalSens)
        results.add_interpolated_pressure_for_physical_sens(interpolatedPressureForMappingPhysicalSens)
        results.add_interpolated_concentration_for_physical_sens(concentrationOfAbsorbanceObservablesForSens)
        return results 

    if physical_sens==1 and bool(observables)and kinetic_sens==1:
        results = simulations.model_data('Shock-Tube',kinetic_sens = S,physical_sens = pS, Solution = timeHistory, Index = ksensIndex, pIndex = psensIndex)
        results.assign_ksens_mappings(Sa ,Sb, Se)
        results.add_interpolated_temps_for_physical_sens(interpolatedTemperatureForMappingPhysicalSens)
        results.add_interpolated_pressure_for_physical_sens(interpolatedPressureForMappingPhysicalSens) 
        results.add_interpolated_concentration_for_physical_sens(concentrationOfAbsorbanceObservablesForSens)
        return results 
#      
    if physical_sens == 0 and kinetic_sens == 0:
        results = simulations.model_data('Shock-Tube',Solution = timeHistory)
        return results