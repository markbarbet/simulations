# -*- coding: utf-8 -*-
"""
Created on Tue Jun 05 15:04:53 2018

@author: HP USER
"""


import numpy as np
import os
import solver as s
import plotter as p

os.chdir('C:\Users\HP USER\Google Drive\Burke Group\Codes')

#########################################################################################
#########################################################################################
#Change inputs here

Temps=np.linspace(800,1100,10)
filenames=[os.getcwd()+'\\Mechanisms\\FFCM-1\\FFCM1.cti','gri30.xml']
pressures=[1]

resTime=[2]
conditions=[{'CH4':1,'O2':1,'N2':3.76}]
observables=['O2','CO']
volume=0.00009
sens=1   #Edit to 1 to run sensitivities
speciesToPlot=['OH','O2']
#########################################################################################
#########################################################################################

results=s.JSR(filenames,Temps,pressures,resTime,conditions,sens,volume,observables,fuel)

a=p.plotSpecies(results,speciesToPlot,conditions)
if sens==1:
    p.plotsens(a,observables,conditions,filenames,Temps,4)