# -*- coding: utf-8 -*-
"""
Created on Fri Feb 02 13:01:11 2018

@author: Mark Barbet
"""
import os
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import troe_rate as tr
gas=ct.Solution(os.getcwd()+'\\TMRP codes\\FFCM-1\\h_c2h2')
kp=[]
kcalc=[]

pressures=np.arange(0.0001*ct.one_atm,100000*ct.one_atm,50000)
pressures=np.arange(0.00001,100000*ct.one_atm,50000)
high_rate=gas.reactions()[0].high_rate.pre_exponential_factor*(1000**gas.reactions()[0].high_rate.temperature_exponent)*np.exp(-gas.reactions()[0].high_rate.activation_energy/1000.0/ct.gas_constant)
low_rate=gas.reactions()[0].low_rate.pre_exponential_factor*(1000**gas.reactions()[0].low_rate.temperature_exponent)*np.exp(gas.reactions()[0].low_rate.activation_energy/1000.0/ct.gas_constant)


for i in np.arange(len(pressures)):
    gas.TP=1000.0,pressures[i]
    kp.append(gas.forward_rate_constants[0])
    kcalc.append(tr.rate(low_rate,high_rate,0.215,10.7,1043.0,2341.0,1000.0,pressures[i]))
f=plt.figure()  
#plt.plot(np.log(pressures),np.log(kp))
plt.loglog(pressures/ct.one_atm,kp,pressures/ct.one_atm,kcalc)
#high_rate=gas.reactions()[0].high_rate.pre_exponential_factor*(1000**gas.reactions()[0].high_rate.temperature_exponent)*np.exp(-gas.reactions()[0].high_rate.activation_energy/1000.0/ct.gas_constant)
#plt.plot(np.log(pressures),np.log(high_rate*np.ones(len(pressures))),'k--') 
plt.loglog(pressures/ct.one_atm,high_rate*np.ones(len(pressures)),'k--')
#low_rate=gas.reactions()[0].low_rate.pre_exponential_factor*(1000**gas.reactions()[0].low_rate.temperature_exponent)*np.exp(gas.reactions()[0].low_rate.activation_energy/1000.0/ct.gas_constant)
temprate=low_rate*pressures/ct.gas_constant/1000.0
#plt.plot(np.log(pressures/ct.one_atm),np.log(rates),'k--')
rates=[]
p2=[]
for i in np.arange(len(temprate)):
    if temprate[i]<high_rate:
        rates.append(temprate[i])
        p2.append(pressures[i])
plt.loglog(np.array(p2)/ct.one_atm,rates,'k--')
plt.xlabel('Pressures (atm)')
plt.ylabel('Rate Constants')
plt.annotate('Temperature = 1000 K',xy=(10e-3,10e1))
plt.savefig(os.getcwd()+'\\figures\\symposiumFigs\\h_c2h2_troe_1000.pdf',dpi=1200,bbox_inches='tight')
low_est=low_rate*ct.one_atm/ct.gas_constant/1000.0
true_est=tr.rate(low_rate,high_rate,0.215,10.7,1043.0,2341.0,1000.0,ct.one_atm,params=['Pr'])
deviation=(low_est-true_est)/true_est


gas=ct.Solution(os.getcwd()+'\\TMRP codes\\FFCM-1\\h_c2h2')
kp=[]
kcalc=[]
pressures=np.arange(0.0001*ct.one_atm,100000*ct.one_atm,50000)
high_rate=gas.reactions()[0].high_rate.pre_exponential_factor*(2000**gas.reactions()[0].high_rate.temperature_exponent)*np.exp(-gas.reactions()[0].high_rate.activation_energy/2000.0/ct.gas_constant)
low_rate=gas.reactions()[0].low_rate.pre_exponential_factor*(2000**gas.reactions()[0].low_rate.temperature_exponent)*np.exp(gas.reactions()[0].low_rate.activation_energy/2000.0/ct.gas_constant)
pressures=np.arange(0.00001,100000*ct.one_atm,50000)


for i in np.arange(len(pressures)):
    gas.TP=2000.0,pressures[i]
    kp.append(gas.forward_rate_constants[0])
    kcalc.append(tr.rate(low_rate,high_rate,0.215,10.7,1043.0,2341.0,2000.0,pressures[i]))

f1=plt.figure()   
#plt.plot(np.log(pressures),np.log(kp))
plt.loglog(pressures/ct.one_atm,kp,pressures/ct.one_atm,kcalc)
#high_rate=gas.reactions()[0].high_rate.pre_exponential_factor*(2000**gas.reactions()[0].high_rate.temperature_exponent)*np.exp(-gas.reactions()[0].high_rate.activation_energy/2000.0/ct.gas_constant)

plt.loglog(pressures/ct.one_atm,high_rate*np.ones(len(pressures)),'k--')
#plt.plot(np.log(pressures),np.log(high_rate*np.ones(len(pressures))),'k--') 
#low_rate=gas.reactions()[0].low_rate.pre_exponential_factor*(2000**gas.reactions()[0].low_rate.temperature_exponent)*np.exp(gas.reactions()[0].low_rate.activation_energy/2000.0/ct.gas_constant)
temprate=low_rate*pressures/ct.gas_constant/2000.0
rates=[]
p2=[]
for i in np.arange(len(temprate)):
    if temprate[i]<high_rate:
        rates.append(temprate[i])
        p2.append(pressures[i])
plt.loglog(np.array(p2)/ct.one_atm,rates,'k--')
#plt.loglog(np.log(pressures/ct.one_atm),np.log(rates),'k--')
plt.xlabel('Pressures (atm)')
plt.ylabel('Rate Constants')
plt.annotate('Temperature = 2000 K',xy=(10e-3,10e1))

plt.savefig(os.getcwd()+'\\figures\\symposiumFigs\\h_c2h2_troe_2000.pdf',dpi=1200,bbox_inches='tight')


