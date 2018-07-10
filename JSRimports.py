# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:23:06 2017

@author: Carly LaGrotta
"""

import yaml
def decisionmakingfunction(filename):
    config = yaml.load(open(filename))
    simulationType = config['apparatus']['kind']
    if simulationType == 'jet stirred reactor' or simulationType=='jsr':
        print(filename)
        a = importingJSR(filename)
        return a
    else: 
        print('we do not have this simulation installed')
        
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
    return(reactorVolume, residenceTime, pressure, initialTemperature,
          finalTemperature, temperatureStep,conditions)