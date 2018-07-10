# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:41:34 2017

@author: Carly LaGrotta
"""

import numpy as np
#do all relevent imports 
import yaml, os
import yaml_parse as dd
import shockTube_constant_pressure as stp
import shockTube_constant_volume as stv
import master_equation as me
import cti_combine as cc
import simulations as sim
import cantera as ct
import free_flame as ff
import experimental_data_import_model_data_object as exp
import SMatrix as sm
   
def simType(filename , absorptionFile = ''):  
    with open(filename) as f:
        config = yaml.load(f)
    simulationType = config['apparatus']['kind']
    experimentType = config['experiment-type']
    mechanicalBoundary = config['common-properties']['assumptions']['mechanical-boundary']

    
    if absorptionFile:
        with open(absorptionFile) as a:
            absorb = yaml.load(a)
        absorption = True
    else:
        absorption = False


            
    possibleSimulationType = ['jet stirred reactor', 'jsr' , 'JSR', 'Jet Stirred Reactor', 
    'shock tube','Shock Tube' ,
    'batch reactor' ,'Batch Reactor',
    'free flame','Free Flame',
    'burner stabalized flame','Burnner Stabalized Flame',
    'ignition delay' ,'Ignition Delay' ,
    'Flame Speeds' ,'flame speeds']
    if simulationType in possibleSimulationType and absorptionFile == '':
        print(simulationType)
        return simulationType,experimentType, config , absorption, mechanicalBoundary
        
    elif simulationType in possibleSimulationType and absorptionFile:
        print(simulationType)
        return simulationType,experimentType,config,absorb,absorption , mechanicalBoundary
        
    else: 
        raise Exception("We do not have this simulation installed, please make sure the name is entered correctly")

        
    
def speci(experimentType):
    if experimentType == 'species profile':
        speciation = True
    
    #how do we want to lable this in the yaml file. Will it be like 'ignition delay'?
    elif experimentType =='ignition delay' or experimentType =='Ignition Delay' or experimentType =='Flame Speeds' or experimentType=='flame speeds': 
        speciation = False
        
    else:
        raise Exception('please check the experiment-type in the yaml file')
        
    return speciation
        
     
    
def Main(filename,ctiFile,absorptionFile='',kineticSens = 0, physicalSens = 0,master_input='',master_index=[]):
    x={}
    if bool(absorptionFile):
        simulationType, experimentType, config, absorb,absorption , mechanicalBoundary = simType(filename , absorptionFile)
    
    else:
        simulationType, experimentType,config,absorption , mechanicalBoundary = simType(filename)
        
    
    
    speciation = speci(experimentType)
    master_cti=me.master_cti(master_input)
    cti_new , originalReactionEquations, masterReactionEquations =cc.cti_write2(x=x,original_cti=ctiFile,master_rxns=master_cti,master_index=master_index)
    print(cti_new)
    #cti_new=ctiFile
    if bool(absorptionFile) and speciation == True and simulationType == 'jsr' or simulationType=='JSR' or simulationType == 'jet stirred reactor' or simulationType == 'Jet Stirred Reactor':
        experimentConditions = dd.importingJSR(config,absorb)
        gas = ct.Solution(cti_new)
        temperatures = experimentConditions['temperature']
       
        modelData = sim.multiTemp(cti_new,gas,temperatures,
                                 experimentConditions,
                                 kinetic_sens = kineticSens,
                                 physical_sens = physicalSens,
                                 observables=['R2OH','H2O','O2'])
                                #physical_params=[])
        print(experimentConditions['observables'])
    elif absorptionFile == '' and speciation == True and simulationType == 'jsr' or simulationType=='JSR' or simulationType == 'jet stirred reactor' or simulationType == 'Jet Stirred Reactor':
        experimentConditions == dd.importingJSR(config)
        gas = ct.Solution(cti_new)
       
       #call multi temp function 
       #gas needs to be defined outsxdie of the 
      
       
        temperatures=experimentConditions['temperature']
        modelData = sim.multiTemp(cti_new,gas,temperatures,
                                 experimentConditions,
                                 kinetic_sens = kineticSens,
                                 physical_sens = physicalSens,
                                 observables=['R2OH','H2O','O2'])
                                #physical_params=[])
       #what are we passing in as physiccal parameters? 
       #where are you getting your species / concentration infromation 
                                     
                                        
                                        
                                        
                                       
       
    elif bool(absorptionFile) and speciation == True and mechanicalBoundary == 'constnat volume' and  simulationType == 'shock tube' or simulationType == 'Shock Tube':
        print(bool(absorptionFile))
        experimentConditions = dd.importingShockTube(config,absorb)
       
        modelData = stv.ShockTube(cti_new, 
                                      experimentConditions['speciesNames'],
                                      experimentConditions['pressure'], 
                                      experimentConditions['temperature'],
                                      experimentConditions['conditions'],
                                      experimentConditions['initialTime'], 
                                      experimentConditions['finalTime'], 
                                      experimentConditions['thermalBoundary'],
                                      experimentConditions['observables'], 
                                      kinetic_sens = kineticSens, 
                                      physical_sens = physicalSens)
        
        if experimentConditions['moleFractionCsvFiles'][0] is not None:
            modelData = exp.importingMoleFractionData(experimentConditions['moleFractionCsvFiles'],experimentConditions['moleFractionObservables'], modelData)
           
        if experimentConditions['absorbanceCsvFiles'][0] is not None:
            modelData = exp.importingAbsorbanceData(experimentConditions['absorbanceCsvFiles'],modelData,absorb, experimentConditions['pathLength'],experimentConditions['absorbanceCsvWavelengths'],kinetic_sens = kineticSens)
            
        if experimentConditions['concentrationCsvFiles'][0] is not None:
             modelData = exp.importingConcentrationData(experimentConditions['concentrationCsvFiles'],experimentConditions['concentrationObservables'], modelData)
        modelData = sm.sMatrixAssembly(modelData,
                                       experimentConditions['moleFractionObservables'],
                                       experimentConditions['concentrationObservables'],
                                       absorbanceObservables = experimentConditions['absorbanceObservables'],
                                       waveLengths = experimentConditions['absorbanceCsvWavelengths'],
                                       kinetic_sens = kineticSens)            
        
        

    elif absorptionFile =='' and speciation == True and mechanicalBoundary == 'constnat volume' and simulationType == 'shock tube' or simulationType == 'Shock Tube':
        experimentConditions = dd.importingShockTube(config)
       
        modelData = stv.ShockTube(cti_new, 
                                      experimentConditions['speciesNames'],
                                      experimentConditions['pressure'], 
                                      experimentConditions['temperature'],
                                      experimentConditions['conditions'],
                                      experimentConditions['initialTime'], 
                                      experimentConditions['finalTime'], 
                                      experimentConditions['thermalBoundary'],
                                      experimentConditions['observables'], 
                                      kinetic_sens = kineticSens, 
                                      physical_sens = physicalSens)  
        if experimentConditions['moleFractionCsvFiles'][0] is not None:
            modelData = exp.importingMoleFractionData(experimentConditions['moleFractionCsvFiles'],experimentConditions['moleFractionObservables'], modelData)
            
        modelData = sm.sMatrixAssembly(modelData,
                                       experimentConditions['moleFractionObservables'],
                                       experimentConditions['concentrationObservables'],
                                       absorbanceObservables = [],
                                       waveLengths = [],
                                       kinetic_sens = kineticSens)            
        if experimentConditions['concentrationCsvFiles'][0] is not None:
             modelData = exp.importingConcentrationData(experimentConditions['concentrationCsvFiles'],experimentConditions['concentrationObservables'], modelData)        
        
    elif absorptionFile =='' and speciation == True and mechanicalBoundary == 'constant pressure' and simulationType == 'shock tube' or simulationType == 'Shock Tube':
        
        experimentConditions = dd.importingShockTube(config)
        modelData = stp.ShockTube(cti_new, 
                                      experimentConditions['speciesNames'],
                                      experimentConditions['pressure'], 
                                      experimentConditions['temperature'],
                                      experimentConditions['conditions'],
                                      experimentConditions['initialTime'], 
                                      experimentConditions['finalTime'], 
                                      experimentConditions['thermalBoundary'],
                                      experimentConditions['observables'], 
                                      kinetic_sens = kineticSens, 
                                      physical_sens = physicalSens)
        if experimentConditions['moleFractionCsvFiles'][0] is not None:
            modelData = exp.importingMoleFractionData(experimentConditions['moleFractionCsvFiles'],experimentConditions['moleFractionObservables'], modelData)
           
            
        if experimentConditions['concentrationCsvFiles'][0] is not None:
             modelData = exp.importingConcentrationData(experimentConditions['concentrationCsvFiles'],experimentConditions['concentrationObservables'], modelData)        
        modelData = sm.sMatrixAssembly(modelData,
                                       experimentConditions['moleFractionObservables'],
                                       experimentConditions['concentrationObservables'],
                                       absorbanceObservables = [],
                                       waveLengths = [],
                                       kinetic_sens = kineticSens)        
        
    elif bool(absorptionFile) and speciation == True and mechanicalBoundary == 'constant pressure' and simulationType == 'shock tube' or simulationType == 'Shock Tube':
        experimentConditions = dd.importingShockTube(config,absorb)

        modelData = stp.ShockTube(cti_new, 
                                      experimentConditions['speciesNames'],
                                      experimentConditions['pressure'], 
                                      experimentConditions['temperature'],
                                      experimentConditions['conditions'],
                                      experimentConditions['initialTime'], 
                                      experimentConditions['finalTime'], 
                                      experimentConditions['thermalBoundary'],
                                      experimentConditions['observables'], 
                                      kinetic_sens = kineticSens, 
                                      physical_sens = physicalSens) 
        
        if experimentConditions['moleFractionCsvFiles'][0] is not None:
            print(experimentConditions['moleFractionCsvFiles'])
            modelData = exp.importingMoleFractionData(experimentConditions['moleFractionCsvFiles'],experimentConditions['moleFractionObservables'], modelData)
            
           
        if experimentConditions['absorbanceCsvFiles'][0] is not None:
            print(experimentConditions['absorbanceCsvFiles'])
            modelData = exp.importingAbsorbanceData(experimentConditions['absorbanceCsvFiles'],modelData,absorb, experimentConditions['pathLength'],experimentConditions['absorbanceCsvWavelengths'],kinetic_sens = kineticSens)
            
        if experimentConditions['concentrationCsvFiles'][0] is not None:
            print(experimentConditions['concentrationCsvFiles'])
            modelData = exp.importingConcentrationData(experimentConditions['concentrationCsvFiles'],experimentConditions['concentrationObservables'], modelData)
        
           
        modelData = sm.sMatrixAssembly(modelData,
                                       experimentConditions['moleFractionObservables'],
                                       experimentConditions['concentrationObservables'],
                                       absorbanceObservables = experimentConditions['absorbanceObservables'],
                                       waveLengths = experimentConditions['absorbanceCsvWavelengths'],
                                       kinetic_sens = kineticSens)
        
        #elif speciation==True and simulationType == 'free flame' or simulationType=='Free Flame' or simulationType=='free-flame' or simulationType=='Free-Flame':
    else:    
        reactants = 'O2:0.23883, N2:0.48466, H2:0.38650, H:0, O:0, OH:0,HO2:0'
        gas=ct.Solution('gri30.xml')
        grid=np.linspace(0,0.02,10)
        mdot = 0.02655058
        gas.TP=400,ct.one_atm
        modelData=ff.free_flame(gas,grid,mdot,kinetic_sens=0,physical_sens=0,energycon=True)
        
    return modelData
#def JSR_isothermal_stdst(Temps,gas,Pressure,concentrations,residenceTime,reactorVolume,pressureValveCoefficient=0.01,maxsimulationTime=1000):
#def JSR_steadystate(gas,resTime,volume,kinetic_sens=0,physical_sens=0,observables=[],physical_params=[],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=1000):
        
    


