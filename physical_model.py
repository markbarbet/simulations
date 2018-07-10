# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 18:41:26 2017

@author: Mark Barbet
"""
import numpy as np
import JSR_funcs as JSR
import JSRimports_new_editing as imp



def run_physical_model(a_params,kinetic_model,filename):
    #YAML reader
    
    if speciation==True:
        run_pmodel_speciation(a_params,kinetic_model,filename)
        
        
        
    else:
        run_pmodel_combustionproperty(a_params,kinetic_model,filename)
        
        
    return model_data_object



def run_pmodel_speciation(a_params,kinetic_model,filename,psens=0,ksens=0):
    
    #Some placeholder stuff for yaml treatment with filename as argument
    
    
    if etype=='JSR' or etype=='jet stirred reactor' or etype =='jet-stirred-reactor':
        f=imp.importingJSR(filename)
        Temps = np.linspace(f['initialTemperature'],f['finalTemperature'],num=f['temperatureStep'])
        gas = ct.Solution(kinetic_model)
        model_data_object=JSR.multiTemp(gas,Temps,f,kinetic_sens=ksens,physical_sens=psens,observables=f['observables'])
        
        
    
    elif etype=='batch reactor' or etype=='batch-reactor':
        
        model_data_object = run_batch_reactor()
        
    elif etype=='free flame' or etype=='free-flame':
        
        model_data_object = run_free_flame(speciation)
        
    elif etype=='burner flame' or etype=='burner-flame':
        
        model_data_object = run_burner_flame(speciation)
        
    elif etype=='shock tube':
        
        model_data_object = shock_tube()
        
    
    
    
    return model_data_object
            
            
def run_pmodel_combustionproperty(a_params,kinetic_model,filename):
    #Placeholder for some yaml reader function with filename as argument
    
    
    if etype=='free flame' or etype=='free-flame':
        model_data_object=run_free_flame_speed()
        
    else if etype=='ignition delay' or etype=='ignition-delay':
        model_data_object=run_ignition_delay()
        
        
        
        
    return model_data_object