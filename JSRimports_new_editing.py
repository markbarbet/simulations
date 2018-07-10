# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:23:06 2017

@author: Carly LaGrotta
"""

import yaml
def sorting(filename):
    with open(filename) as f:
        config = yaml.load(f)
    simulationType = config['apparatus']['kind']
    if simulationType == 'jet stirred reactor' or simulationType=='jsr':
        print(simulationType)
        a = importingJSR(filename)
        return a
    else: 
        raise Exception("We do not have this simulation installed")
        
def importingJSR(filename):
    with open(filename) as f:
        config = yaml.load(f)
    
    reactorVolume = config['apparatus']['reactor-volume']['value']
    residenceTime = config['apparatus']['residence-time']['value']
    pressure = config['common-properties']['pressure']['value']
    initialTemperature = config['common-properties']['temperature']['initial-value']
    finalTemperature = config['common-properties']['temperature']['final-value']
    temperatureStep = config['common-properties']['temperature']['step']
    moleFractions = [((concentration['mole-fraction'])) for concentration in config['common-properties']['composition']]
    speciesNames = [(species['species']) for species in config['common-properties']['composition']]
    conditions = dict(zip(moleFractions,speciesNames))
    moleFractionObservables = [datapoint['targets'][0]['name'] for datapoint in config['datapoints']['mole-fraction']]
    absorbanceObservables = [species['species'] for species in config['datapoints']['absorbance']['absorbing-species']]
    observables = moleFractionObservables + absorbanceObservables
    return {
            'reactorVolume': reactorVolume,
            'residenceTime': residenceTime,
            'pressure': pressure,
            'initialTemperature': initialTemperature,
            'finalTemperature': finalTemperature,
            'temperatureStep': temperatureStep,
            'conditions': conditions,
            'observables': observables
    }

def importingShockTube(filename):
    with open(filename) as f:
        config = yaml.load(f)
    pressure = config['common-properties']['pressure']['value']
    temperature = config['common-properties']['temperature']['value']
    moleFractions = [((concentration['mole-fraction'])) for concentration in config['common-properties']['composition']]
    speciesNames = [(species['species']) for species in config['common-properties']['composition']]
    conditions = dict(zip(moleFractions,speciesNames))
    
    return{
           'pressure' : pressure,
           'temperature' : temperature,
           'conditions' : conditions
           
              
           }
    
  