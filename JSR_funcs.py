# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:57:52 2017

@author: HP USER
"""
from __future__ import division
from __future__ import print_function
import numpy as np
import cantera as ct
import time
import pandas as pd


class jsr_model_data:
    def __init__(self,kinetic_sens=[],physical_sens=[],reactorSolution=[]):
        self.k_sens=kinetic_sens
        self.p_sens=physical_sens
        self.speciation=reactorSolution





def JSR_isothermal_stdst(Temps,gas,Pressure,concentrations,residenceTime,reactorVolume,pressureValveCoefficient=0.01,maxsimulationTime=1000):
    stirredReactor_list = []
    gas_List=[]
    for i in np.arange(len(Temps)):
        # Inlet gas conditions
        reactorTemperature = Temps[i] #Kelvin
        reactorPressure = Pressure #in atm.
        gas.TPX = reactorTemperature, reactorPressure, concentrations
        # Reactor parameters
        # Instrument parameters

        # This is the "conductance" of the pressure valve and will determine its efficiency in 
        # holding the reactor pressure to the desired conditions. 
        pressureValveCoefficient=pressureValveCoefficient

        # This parameter will allow you to decide if the valve's conductance is acceptable. If there
        # is a pressure rise in the reactor beyond this tolerance, you will get a warning
        maxPressureRiseAllowed = 0.001

        # Simulation termination criterion
        maxSimulationTime = maxsimulationTime # seconds
        fuelAirMixtureTank = ct.Reservoir(gas)
        exhaust = ct.Reservoir(gas)

        stirredReactor = ct.IdealGasReactor(gas, energy='off', volume=reactorVolume)

        massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,downstream=stirredReactor,mdot=stirredReactor.mass/residenceTime)

        pressureRegulator = ct.Valve(upstream=stirredReactor,downstream=exhaust,K=pressureValveCoefficient)

        reactorNetwork = ct.ReactorNet([stirredReactor])

        # now compile a list of all variables for which we will store data
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        columnNames = ['pressure'] + columnNames

        # use the above list to create a DataFrame
        timeHistory = pd.DataFrame(columns=columnNames)

        # Start the stopwatch
        tic = time.time()

        # Set simulation start time to zero
        t = 0
        counter = 1
        while t < maxSimulationTime:
            t = reactorNetwork.step()

            # We will store only every 10th value. Remember, we have 1200+ species, so there will be
            # 1200 columns for us to work with
            if(counter%10 == 0):
            #Extract the state of the reactor
                state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                   stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
        
            #Update the dataframe
                timeHistory.loc[t] = state
    
            counter += 1

        # Stop the stopwatch
        toc = time.time()

        print('Simulation Took {:3.2f}s to compute, with {} steps'.format(toc-tic, counter))
        #state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                   #stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
        #for j in np.arange(len(speciesName)):
            #component_final_X[i,j] = stirredReactor.get_state()[stirredReactor.component_index(speciesName[j])]
            #component_final_X[i,j] = state[stirredReactor.component_index(speciesName[j])+1]
        # We now check to see if the pressure rise during the simulation, a.k.a the pressure valve
        # was okay
        pressureDifferential = timeHistory['pressure'].max()-timeHistory['pressure'].min()
        if(abs(pressureDifferential/reactorPressure) > maxPressureRiseAllowed):
            print("WARNING: Non-trivial pressure rise in the reactor. Adjust K value in valve")
    
        stirredReactor_list.append(stirredReactor)
        
    #%matplotlib notebook

    #plt.style.use('ggplot')
    #plt.style.use('seaborn-pastel')

    #plt.rcParams['axes.labelsize'] = 18
    #plt.rcParams['xtick.labelsize'] = 14
    #plt.rcParams['ytick.labelsize'] = 14
    #plt.rcParams['figure.autolayout'] = True

    #plt.figure()
    #plt.semilogx(timeHistory.index, timeHistory['B2CO'],'-o')
    #plt.xlabel('Time (s)')
    #plt.ylabel(r'Mole Fraction : $X_{CO}$');
    return stirredReactor_list, gas

def JSR_steadystate(gas,resTime,volume,kinetic_sens=0,physical_sens=0,observables=[],physical_params=[],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=1000):
    # Inlet gas conditions are passed into function in the "gas" object, which is a cantera object
    # Reactor parameters passed into function as resTime and volume.  Residence time and volume of JSR
    #kinetic sens and physical sens are optional parameters which are set to zero by default.  Set them to 1 to
    #calculate sensitivity based on kinetic or physical parameters.  If these are set to 1 you must pass
    #an array of all observables to calculate sensitivities for
    reactorPressure=gas.P
    # This is the "conductance" of the pressure valve and will determine its efficiency in 
    # holding the reactor pressure to the desired conditions. It is an optional parameter
    pressureValveCoefficient=pressureValveCoefficient

    # This parameter will allow you to decide if the valve's conductance is acceptable. If there
    # is a pressure rise in the reactor beyond this tolerance, you will get a warning
    maxPressureRiseAllowed = 0.001

    # Simulation termination criterion
    maxSimulationTime = maxsimulationTime # seconds.  An optional parameter
    fuelAirMixtureTank = ct.Reservoir(gas)
    exhaust = ct.Reservoir(gas)

    stirredReactor = ct.IdealGasReactor(gas, energy=energycon, volume=volume)
    
    
        
    massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,downstream=stirredReactor,mdot=stirredReactor.mass/resTime)

    pressureRegulator = ct.Valve(upstream=stirredReactor,downstream=exhaust,K=pressureValveCoefficient)

    reactorNetwork = ct.ReactorNet([stirredReactor])
    
    #This block adds kinetic sensitivity parameters for all reactions if desired.  
    if kinetic_sens==1 and bool(observables):
        for i in range(gas.n_reactions):
            stirredReactor.add_sensitivity_reaction(i)
    if kinetic_sens==1 and bool(observables)==False:
        print('Please supply a non-empty list of observables for sensitivity analysis or set kinetic_sens=0')
        
    if physical_sens==1 and bool(observables):
        print('Placeholder')
    if physical_sens==1 and bool(observables)==False:
        print('Please supply a non-empty list of observables for sensitivity analysis or set physical_sens=0')

    # now compile a list of all variables for which we will store data
    columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
    columnNames = ['pressure'] + columnNames

    # use the above list to create a DataFrame
    timeHistory = pd.DataFrame(columns=columnNames)

    # Start the stopwatch
    tic = time.time()
        
    #Establish a matrix to hold sensitivities for kinetic parameters, along with tolerances
    if kinetic_sens==1 and bool(observables):
        senscolumnNames = ['Reaction']+observables      
        sensArray = pd.DataFrame(columns=senscolumnNames)
        senstempArray = np.zeros((len(observables),gas.n_reactions))
        reactorNetwork.rtol_sensitivity = 1.0e-6
        reactorNetwork.atol_sensitivity = 1.0e-6
    if physical_sens==1 and bool(observables):
        senscolumnNames = ['Parameter'] + observables
        psensArray = pd.DataFrame(columns=senscolumnNames)
        senstempArray = np.zeros((len(observables),len(physical_params)))
        
    reactorNetwork.advance_to_steady_state()
    
    if kinetic_sens==1 and bool(observables):
        for k in np.arange(len(observables)):
            for j in np.arange(gas.n_reactions):
                    #print( k,j)
                senstempArray[k,j]=reactorNetwork.sensitivity(observables[k],j)
        sensArray['Reaction']=gas.reaction_equations()        
        sensArray[observables]=senstempArray.T   
    
    
    #state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                   #stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X]) 
        
        # Stop the stopwatch
    toc = time.time()

    print('Simulation Took {:3.2f}s to compute'.format(toc-tic))
    
    # We now check to see if the pressure rise during the simulation, a.k.a the pressure valve
    # was okay
    pressureDifferential = timeHistory['pressure'].max()-timeHistory['pressure'].min()
    if(abs(pressureDifferential/reactorPressure) > maxPressureRiseAllowed):
        print("WARNING: Non-trivial pressure rise in the reactor. Adjust K value in valve")
        
    if kinetic_sens==1 and bool(observables) and physical_sens==0:
        modelData=jsr_model_data(kinetic_sens=sensArray,reactorSolution=stirredReactor)
        return modelData
    if physical_sens==1 and bool(observables)and kinetic_sens==0:
        modelData=jsr_model_data(physical_sens=psensArray,reactorSolution=stirredReactor)
    if physical_sens==1 and bool(observables) and kinetic_sens==1:
        modelData=jsr_model_data(physical_sens=psensArray,kinetic_sens=sensArray,reactorSolution=stirredReactor)
    else:
        modelData=jsr_model_data(reactorSolution=stirredReactor)
        return modelData
    
    
def multiTemp(gas,Temps,f,kinetic_sens=0,physical_sens=0,observables=[],physical_params=[],energycon='off',pressureValveCoefficient=0.01,maxsimulationTime=1000):
    modelDataList = [] 
    for i in np.arange(len(Temps)):
        gas.TPX = Temps[i],f['pressure'],f['conditions']
        solutionObject = JSR_steadystate(gas,f['residenceTime'],f['reactorVolume'])
        modelDataList.append(solutionObject)
    