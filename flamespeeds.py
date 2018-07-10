# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 11:53:46 2017

@author: Mark Barbet
"""

#Free Flame simulator for use in TMRP project

import numpy as np
import cantera as ct
import pandas as pd
import os
import free_flame as ff
#import simulations as sim

fuels=['CH3OH']
#fuels=[{'H2':1,'CO':1}]
phi=[]
phi=[1]
phi=[0.6,0.9,1.0,1.2,1.4,1.6,1.8]
phi=np.arange(0.6,2.0,0.2)
#phi=[0.4]
#phi=[0.6]
#results=[]
if 'results' not in locals().keys():
    results=[]

pressures=[1]
temps=[298]

#mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta1_forward_specified2500.0_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\JetSurf\\jetsurf2.cti',
            os.getcwd()+'\\TMRP codes\\JetSurf\\appendedinput_beta001_forward_chem.cti']

#mechanisms=[os.getcwd()+'\\TMRP codes\\dooley\\chem.cti',
#            os.getcwd()+'\\TMRP codes\\dooley\\appendedinput_beta001_forward_chem.cti',
#            os.getcwd()+'\\TMRP codes\\dooley\\appendedinput_beta001_chem.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_FFCM1.cti',            
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_specified5000.0_top30_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_top30_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top30_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_low_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_low_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_C2H2_forward_focus_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_cutrxns01_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_cutrxns01_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_top29_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top29_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_top30_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top30_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1.cti',
            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_minusCutrxns_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_modify948_FFCM1.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\Metcalfe\\chem.cti',
            os.getcwd()+'\\TMRP codes\\Metcalfe\\appendedinput_beta001_forward_top30_chem.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_FFCM1_H.cti',            
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_specified5000.0_FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_specified5000.0_top6_FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_top6_FFCM1_H.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_top6_FFCM1_H.cti']
mechanisms=[os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem.cti']
mechanisms=[os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1.cti',
            os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1_natureTermoleculars.cti',
            os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1_natureTermolecularsx10.cti']
width=0.3

for i in np.arange(len(phi)):
    for filename in np.arange(len(mechanisms)):
        
        oxidizer={'O2':0.5, 'N2':0.5*3.76}
        
        gas=ct.Solution(mechanisms[filename])
        gas.TP=temps[0],pressures[0]*ct.one_atm
        results.append(ff.free_flame(phi[i],fuels[0],oxidizer,gas,width,kinetic_sens=0,energycon=True,flamespeed_sens=1,soret=False))
        results[-1].add_mechanism(mechanisms[filename])
#results[0].plot_flamespeed_sens(1,'appendedinput_FFCM1',0.01)
#runfile('C:/Users/HP USER/Google Drive/Burke Group/Codes/TMRP_analysis.py', wdir='C:/Users/HP USER/Google Drive/Burke Group/Codes')
        
import shelve
import os
#T='Hiya'
#val=[1,2,3]

filename=os.getcwd()+'\\tmp\\flamespeeds_'+results[0].fuel+'methanolTest.out'
my_shelf = shelve.open(filename,'n') # 'n' for new

for key in dir():
    if key=='results':
        print(key)
        my_shelf[key] = globals()[key]
    #except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        #print('ERROR shelving: {0}'.format(key))
my_shelf.close()
#To restore:
