# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 13:27:55 2017

@author: Mark Barbet
"""


import numpy as np
import cantera as ct
import pandas as pd
import os
#import free_flame as ff
import simulations as sim
fuel='H'
fuel='C4H6'
#fuel='CH3OCHO'
phi=1.8
#phi=1
conditions=[{fuel:35.8,'O2':35.8,'AR':24.6}]
#conditions=[{fuel:0.076*phi,'O2':0.231,'AR':0.693}]
#conditions=[{fuel:phi,'O2':0.5,'N2':1.88,'H':0, 'O':0, 'OH':0,'HO2':0}]
#conditions=[{'CH4':1.2,'O2':2,'N2':7.52}]
fuel='H2'
phi=[1.5]
phi=1.5
results=[]
pressures=[0.102632]
pressures=[0.05]
temps=[400]
grid=0.018
sens=0
conditions=[{fuel:phi*1,'O2':0.5,'N2':0.5*3.76,'H':0.0,'O':0.0,'OH':0.0,'HO2':0.0}]
#conditions=[{fuel:7.6,'O2':23.1,'AR':69.3}]
conditions = [{fuel:phi*1,'O2':0.5,'N2':0.5*3.76}]
obs=['O2','CO','CO2','H','H2','C4H6','OH','H2O']
obs=['CO','CO2']
obs=['C4H6','C6H6']
#obs=['CO','CO2','H2O','O2','H2']
#obs=['CO2','CO','O2','CH3OCHO','H2','H2O']
obs=['NO','OH']
velocity=[0.693]
Area=0.25*np.pi*(0.06**2)
Q=0.0001
#Q=5.066667e-5
flowRatebool=True  #This value is true if mdot is calculated from a volumetric flowrate
datafiles=[os.getcwd()+'\\paperFigs\\harrington_38torr_T_profile_dense.csv']
#mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_FFCM1.cti',            
#            os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta001_forward_FFCM1.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\FFCM-1\\appendedinput_beta1_forward_specified2500.0_FFCM1.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti',
##            os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta001_forward_uscmech.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti',
#            os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta001_forward_uscmech.cti',
#            os.getcwd()+'\\TMRP codes\\USCMech\\appendedinput_beta001_uscmech.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\dooley\\chem.cti',
##            os.getcwd()+'\\TMRP codes\\dooley\\appendedinput_beta001_forward_chem.cti',
##            os.getcwd()+'\\TMRP codes\\dooley\\appendedinput_beta001_chem.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\dooley\\chem.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\Lamoreaux\\lam_onlyNO.cti',
##            os.getcwd()+'\\TMRP codes\\Lamoreaux\\appendedinput_beta1_lam_onlyNO.cti',
##            os.getcwd()+'\\TMRP codes\\Lamoreaux\\appendedinput_beta1_forward_lam_onlyNO.cti',
##            os.getcwd()+'\\TMRP codes\\Lamoreaux\\lam_onlyNO_plusTer_beta0.1.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\Lamoreaux\\lam_onlyNO.cti',
##            os.getcwd()+'\\TMRP codes\\Lamoreaux\\lam_onlyNO_plusTer_beta0.1.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\nominal_78Torr\\chem_onlyNO.cti',
##            os.getcwd()+'\\TMRP codes\\nominal_78Torr\\chem_onlyNO_plusTer.cti']
##mechanisms=['gri30.xml']
##mechanisms=['C:\\Users\\HP USER\\Google Drive\\Burke Group\\Codes\\TMRP codes\\konnov0.6\\chem.cti']
#mechanisms=[os.getcwd()+'\\TMRP codes\\konnov0.6\\chem_onlyNO_klipp.cti',
#            os.getcwd()+'\\TMRP codes\\konnov0.6\\appendedinput_beta3_chem_onlyNO_klipp.cti',
#            os.getcwd()+'\\TMRP codes\\konnov0.6\\appendedinput_beta3_forward_chem_onlyNO_klipp.cti']
##mechanisms=[os.getcwd()+'\\TMRP codes\\konnov0.6\\chem_onlyNO_klipp.cti',
#  #          os.getcwd()+'\\TMRP codes\\konnov0.6\\appendedinput_N2O_H_chem_onlyNO.cti']
mechanisms=[os.getcwd()+'\\Mechanisms\\konnov0.6\\chem_onlyNO.cti',
            os.getcwd()+'\\Mechanisms\\Lamoreaux\\lam_onlyNO.cti',
            os.getcwd()+'\\Mechanisms\\glarborg\\chem_onlyNO.cti',
            os.getcwd()+'\\Mechanisms\\glarborg\\chem_deleted.cti']
mechanisms=[os.getcwd()+'\\Mechanisms\\glarborg\\chem_deleted.cti']
with open(datafiles[0]) as f:
    profile=f.readlines()
    
for i in np.arange(len(profile)):
    profile[i]=profile[i].rstrip('\n').rstrip(',')+'\n'

with open(datafiles[0],'w') as f:
    f.writelines(profile)
    
data=pd.read_csv(datafiles[0],header=3)

for i in np.arange(len(conditions)):
    for filename in np.arange(len(mechanisms)):
        
        #oxidizer={'O2':15.5, 'N2':58.28}
        tempgas=ct.Solution(mechanisms[filename])
        tempgas.TPX=300,pressures[0]*ct.one_atm,conditions[i]
        tempgas.TPX=300,ct.one_atm,conditions[i]
        density=tempgas.DP[0]
        
        mdot=density*velocity[0]
        #mdot=0.03047*density
        print(mdot)
        #grid=np.linspace(0,0.03,10)
        #grid=0.02
        if flowRatebool:
            
            m=Q*tempgas.density
            mdot=m/Area
            #mdot=0.026
            print(mdot)
        
        gas=ct.Solution(mechanisms[filename])
        gas.TPX=temps[0],pressures[0]*ct.one_atm,conditions[i]
        #results.append(ff.free_flame(phi[i],fuels,oxidizer,gas,width,kinetic_sens=0,energycon=True,flamespeed_sens=1,soret=False))
        results.append(sim.burner_flame(gas,grid,mdot,data=data,kinetic_sens=sens,observables=obs,soret=True))
        results[-1].add_fuel(fuel)
        results[-1].add_mechanism(mechanisms[filename])
        results[-1].assign_phi(phi)
#results[0].plot_flamespeed_sens(1,'appendedinput_FFCM1',0.01)

#runfile('C:/Users/HP USER/Google Drive/Burke Group/Codes/burnerflame_analysis.py', wdir='C:/Users/HP USER/Google Drive/Burke Group/Codes')