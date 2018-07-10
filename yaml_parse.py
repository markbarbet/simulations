# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:23:06 2017

@author: Carly LaGrotta
"""

import yaml
import numpy as np
def sorting(filename):
    with open(filename) as f:
        config = yaml.load(f)
    simulationType = config['apparatus']['kind']
    if simulationType == 'jet stirred reactor' or simulationType=='jsr' or simulationType=='Jet Stirred Reactor' or simulationType=='shock tube' or simulationType=='Shock Tube' or simulationType=='batch reactor' or simulationType=='Batch Reactor' or simulationType=='free flame'or simulationType=='Free Flame' or  simulationType=='burner flame' or simulationType=='Burnner Flame'or simulationType=='ignition delay' or simulationType=='Ignition Delay':
        print(simulationType)
        return simulationType
    else: 
        raise Exception("We do not have this simulation installed, please make sure the name is entered correctly")
        
def importingJSR(config , absorb = {}):
#    with open(filename) as f:
#        config = yaml.load(f)
    if absorb == {}:    
        reactorVolume = config['apparatus']['reactor-volume']['value']
        residenceTime = config['apparatus']['residence-time']['value']
        pressure = config['common-properties']['pressure']['value']
        temperatures=config['common-properties']['temperature']['values']
        moleFractions = [((concentration['mole-fraction'])) for concentration in config['common-properties']['composition']]
        speciesNames = [(species['species']) for species in config['common-properties']['composition']]
        conditions = dict(zip(speciesNames,moleFractions))
        print(conditions)
        moleFractionObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['mole-fraction']]
        concentrationObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['concentration']]
    #    absorbanceObservables = []
    #    for test in xrange(len(config['datapoints']['absorbance'])):
    #        temp = [spec['species'] for spec in config['datapoints']['absorbance'][test]['absorbing-species']]
    #        absorbanceObservables = absorbanceObservables + temp 
        observables = [x for x in (moleFractionObservables + concentrationObservables) if x is not None]
        return {
                'reactorVolume': reactorVolume,
                'residenceTime': residenceTime,
                'pressure': pressure,
                'temperature':temperatures, 
                'conditions': conditions,
                'observables': observables
                }
    elif absorb:
        reactorVolume = config['apparatus']['reactor-volume']['value']
        residenceTime = config['apparatus']['residence-time']['value']
        pressure = config['common-properties']['pressure']['value']
        temperatures=config['common-properties']['temperature']['values']
        moleFractions = [((concentration['mole-fraction'])) for concentration in config['common-properties']['composition']]
        speciesNames = [(species['species']) for species in config['common-properties']['composition']]
        conditions = dict(zip(speciesNames,moleFractions))
        print(conditions)
        moleFractionObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['mole-fraction']]
        concentrationObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['concentration']]
        
        absorptionObservables = [species['species'] for species in absorb['Absorption-coefficients']]
        observables = [x for x in (moleFractionObservables + concentrationObservables + absorptionObservables) if x is not None]

        
        return {'reactorVolume': reactorVolume,
                'residenceTime': residenceTime,
                'pressure': pressure,
                'temperature':temperatures, 
                'conditions': conditions,
                'observables': observables
                }
        

       
        
    

        
            
        

def importingShockTube(config, absorb = {}):
#    with open(filename) as f:
#        config = yaml.load(f)
    if absorb == {}:
        pressure = config['common-properties']['pressure']['value']
        temperature = config['common-properties']['temperature']['value']
        moleFractions = [((concentration['mole-fraction'])) for concentration in config['common-properties']['composition']]
        speciesNames = [(species['species']) for species in config['common-properties']['composition']]
        conditions = dict(zip(speciesNames,moleFractions))
        thermalBoundary = config['common-properties']['assumptions']['thermal-boundary']
        moleFractionObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['mole-fraction']]

 

        concentrationObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['concentration']]            
        observables = [x for x in (moleFractionObservables + concentrationObservables) if x is not None]

        initialTime = config['common-properties']['time']['initial-time']['value']
        #eventually going to get this from a csv file 
        finalTime = config['common-properties']['time']['final-time']['value']
        moleFractionCsvFiles = [csvfile['csvfile'] for csvfile in config['datapoints']['mole-fraction']]
        absorbanceCsvFiles = [csvfile['csvfile'] for csvfile in config['datapoints']['absorbance']]
        concentrationCsvFiles = [csvfile['csvfile'] for csvfile in config['datapoints']['concentration']]
        csvFiles = [x for x in (moleFractionCsvFiles + absorbanceCsvFiles + concentrationCsvFiles) if x is not None]
    
    
        
        return{
               'pressure' : pressure,
               'temperature' : temperature,
               'conditions' : conditions,
               'thermalBoundary' : thermalBoundary,
               'observables': observables,
               'initialTime': initialTime,
               'finalTime' : finalTime,
               'speciesNames': speciesNames,
               
               
                  
               }
            
    elif absorb:
        
        
        pressure = config['common-properties']['pressure']['value']
        temperature = config['common-properties']['temperature']['value']
        moleFractions = [((concentration['mole-fraction'])) for concentration in config['common-properties']['composition']]
        speciesNames = [(species['species']) for species in config['common-properties']['composition']]
        conditions = dict(zip(speciesNames,moleFractions))
        thermalBoundary = config['common-properties']['assumptions']['thermal-boundary']
        moleFractionObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['mole-fraction']]

#        absorbanceObservables = []
#        absorbanceSpeciesLst = [[] for i in xrange(len(config['datapoints']['absorbance']))]
#        for test in xrange(len(config['datapoints']['absorbance'])):
#            temp = [spec['species'] for spec in config['datapoints']['absorbance'][test]['absorbing-species']]
#            absorbanceSpeciesLst[test] = temp        
#            absorbanceObservables = absorbanceObservables + temp 
        absorptionObservables = [species['species'] for species in absorb['Absorption-coefficients']]
        

        concentrationObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['concentration']]           
        observables = [x for x in (moleFractionObservables + concentrationObservables + absorptionObservables) if x is not None]
        initialTime = config['common-properties']['time']['initial-time']['value']
        #eventually going to get this from a csv file 
        finalTime = config['common-properties']['time']['final-time']['value']
        moleFractionCsvFiles = [csvfile['csvfile'] for csvfile in config['datapoints']['mole-fraction']]
        absorbanceCsvFiles = [csvfile['csvfile'] for csvfile in config['datapoints']['absorbance']]
        concentrationCsvFiles = [csvfile['csvfile'] for csvfile in config['datapoints']['concentration']]
        csvFiles = [x for x in (moleFractionCsvFiles + absorbanceCsvFiles + concentrationCsvFiles) if x is not None]
        
        return {
               'pressure' : pressure,
               'temperature' : temperature,
               'conditions' : conditions,
               'thermalBoundary' : thermalBoundary,
               'observables': observables,
               'initialTime': initialTime,
               'finalTime' : finalTime,
               'speciesNames': speciesNames
               
                  
               } 
        
  