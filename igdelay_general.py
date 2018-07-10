# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:04:47 2018

@author: Mark Barbet
"""

import cantera as ct
import numpy as np
import os
import efficiency_manipulate as em
import ig_delay as ig
import matplotlib.pyplot as plt


#Range of temperatures to run simulations for
T=np.arange(1000,1250,50)#methanol
#T=np.arange(625,1450,100)#dme
val='all'
print(val)
#Pressure input *******IN ATMOSPHERES ONLY*******
P=19.3#methanol
#P=11.0#dme
X={'CH3OH':5.0, 'O2':15.0, 'N2':80}#methanol
#X={'ch3och3':0.3,'o2':3.0,'n2':3.76*3.0} #dme
#P=2.36#propane
#T=np.arange(1300,1500,40)#propane
#X={'C3H8':1.6,'O2':8.0,'AR':90.4}#propane
mechanism=[os.getcwd()+'\\TMRP codes\\dryer_methanol\\chem.cti',os.getcwd()+'\\TMRP codes\\FFCM-1\\FFCM1.cti',os.getcwd()+'\\TMRP codes\\USCMech\\uscmech.cti']
mechanism=[os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem.cti']
#mechanism=[os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem.cti',os.getcwd()+'\\Mechanisms\\dryer_methanol\\chem_methanolEff.cti']

#mechanism=[os.getcwd()+'\\Mechanisms\\CH4_DME\\chem.cti']
mechanism=[os.getcwd()+'\\Mechanisms\\Aramco2.0\\chem.cti']
mechanism=[os.getcwd()+'\\Mechanisms\\nheptane\\nheptane.cti']
mechanism=[os.getcwd()+'\\Mechanisms\\cyclohexane\\chem.cti']
mechanism=[os.getcwd()+'\\Mechanisms\\isopentanol\\chem.cti']
mechanism=[os.getcwd()+'\\Mechanisms\\UdatedH2Model110725\\chem.cti']

P=[1.0,10.0,100.0]
phi=[0.5,1.0,2.0,4.0]
T=np.arange(500.0,2000.0,100)
Fuels=['CH4','C2H6','C2H4','C2H2','C3H8','C3H6','C3H4','C4H10','IC4H10','C4H8-1','C4H8-2','IC4H8','C4H6-1','C4H6-2','C4H6','C4H612']
Fuels=['CH4','C2H6','C2H4']
#Fuels=['NC7H16']
#Fuels=['chx']
#Fuels=['ic5h11oh']
#Fuels=['H2']
manipulated=False
num=[2.0,3.5,3.0,2.5,5.0,4.5,4.0,6.5,6.5,6.0,6.0,6.0,5.5,5.5,5.5,5.5]
num=[2.0,3.5,3.0]
#num=[11.0]
#num=[9.0]
#num=[7.5]
#num=[1.0]
for p in P:
    for x in np.arange(len(phi)):
        for f in np.arange(len(Fuels)):
            files=[]
            files=mechanism
            
            X={Fuels[f]:phi[x],'O2':num[f],'N2':3.76*num[f]}
            efficiency_manipulate=True
            reactorType='cv'
            
            if efficiency_manipulate and manipulated==False:
                gases=[]
                for i in np.arange(len(files)):
                    gas=ct.Solution(files[i])
                    gases.append(gas)
                    gas.name='igDelayRateSwap_'+files[i].split('\\')[-1].rstrip('.cti')+Fuels[f]
                    gas2=em.efficiency_rate_swap(gas,[val])
                    gases.append(gas2)
                    import soln2cti as ctiw
                    new_file=ctiw.write(gas)
                    #gas2=ct.Solution(os.getcwd()+'\\pym_gas.cti')
                    files.append(new_file)
                    manipulated=True
            
            delays=[]
            for i in np.arange(len(files)):
                #gas=ct.Solution(mechanism[i])
                delays.append([])
                for j in np.arange(len(T)):
                    try:
                        delays[i].append(ig.ignition_delay(files[i],T[j],p,X,options=reactorType))
                    except:
                        delays[i].append('Temperature '+str(T[j])+'K failed')
            plt.figure()
            colors=['b','r','k','g','c','m']
            for i in np.arange(len(files)):
                tempdel=[]
                tempT=[]
                for j in np.arange(len(delays[i])):
                    if str(type(delays[i][j]))!="<type 'str'>":
                        tempdel.append(delays[i][j])
                        tempT.append(T[j])
                plt.semilogy(1000.0/np.array(tempT),tempdel,color=colors[i])
                plt.title(Fuels[f]+' at phi= '+str(phi[x])+', '+str(p)+' atm')
    

    


