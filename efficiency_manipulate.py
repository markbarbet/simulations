# -*- coding: utf-8 -*-
"""
Created on Thu Jan 04 19:33:28 2018

@author: Mark Barbet
"""

#Code here is written to perform some changes to the efficiencies in a cantera model to see what differences are caused by how efficiencies 
#are generally defined in most three body and pressure dependent reactions\


import cantera as ct
import numpy as np


#Call this function with a cantera gas object
def efficiency_rate_swap(gas,rxn_indices=[],non_ignored_species=['N2','AR','HE','n2','ar','he','Ar','He','o2','O2']):

    #Find reactions in mechanism that have a collider 'M' that requires efficiency definitions
    if rxn_indices==['all']:
        M_indices=[]
    
        for i in np.arange(len(gas.reactions())):
            if isinstance(gas.reactions()[i],ct._cantera.ThreeBodyReaction) or isinstance(gas.reactions()[i],ct._cantera.FalloffReaction):
                M_indices.append(i)
    elif rxn_indices!=['all'] and rxn_indices!=['none']:
        M_indices=rxn_indices
    elif rxn_indices==['none']:
        M_indices=[]
    #print(M_indices)
    
    if M_indices!=[]:
        for i in M_indices:
            #Find the max efficiency specified for the reaction
            #print(i)
            try:
                for e in np.arange(len(non_ignored_species)):
                    if non_ignored_species[e] in gas.species_names and non_ignored_species[e] not in gas.reactions()[i].efficiencies.keys():
                        
                        temp_dict=gas.reactions()[i].efficiencies
                        temp_dict[non_ignored_species[e]]=1.0
                        gas.reactions()[i].efficiencies=temp_dict
                        temp_dict={}
                
                max_eff = np.max(gas.reactions()[i].efficiencies.values())
                
                #Now write the efficiencies as the original efficiency divided by the max
                #and then multiply the reaction rate constant by the max efficiency
                tempDict=gas.reactions()[i].efficiencies
                for j in gas.reactions()[i].efficiencies:
                    tempDict[j] = np.divide(gas.reactions()[i].efficiencies[j],max_eff)
                    gas.reactions()[i].efficiencies=tempDict
                    
                if isinstance(gas.reactions()[i],ct._cantera.ThreeBodyReaction):
                    gas.reactions()[i].rate=ct.Arrhenius(gas.reactions()[i].rate.pre_exponential_factor*max_eff,gas.reactions()[i].rate.temperature_exponent,
                                  gas.reactions()[i].rate.activation_energy)
                    #print('true')
                elif isinstance(gas.reactions()[i],ct._cantera.FalloffReaction):
                    gas.reactions()[i].high_rate=ct.Arrhenius(gas.reactions()[i].high_rate.pre_exponential_factor*max_eff,gas.reactions()[i].high_rate.temperature_exponent,
                                  gas.reactions()[i].high_rate.activation_energy)
                    gas.reactions()[i].low_rate=ct.Arrhenius(gas.reactions()[i].low_rate.pre_exponential_factor*max_eff,gas.reactions()[i].low_rate.temperature_exponent,
                                  gas.reactions()[i].low_rate.activation_energy)
            except:
                print('Reaction '+str(i)+' may not have efficiencies to manipulate')
    return gas
#Now, the user may simply use this function to operate on any gas solution object of their choice