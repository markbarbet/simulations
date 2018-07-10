# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:53:36 2018

@author: Mark Barbet
"""

#Function to find ignition delays for adiabatic, constant pressure


import cantera as ct
import numpy as np
import os
import shockTube_constant_volume as stcv
import shockTube_constant_pressure as stcp

def ignition_delay(mechanism,T,P,X,options='cv',igdelay_type='T'):
    if options=='cp':
        tt = []
        TT = []
        results=stcp.ShockTube(mechanism,['N2'],P,T,X,0,0.02,'adiabatic')
        #t = 0.0
        #while t < 0.02:
            #t = sim.step()
            #tt.append(1000 * t)
            #TT.append(r.T)
        delay=[]
        tt=results.solution['time'].values
        TT=results.solution['temperature'].values
        dTdt=np.zeros(np.array(tt).shape,np.float)
        dTdt[0:-1]=np.diff(TT)/np.diff(tt)
        dTdt[-1]=(TT[-1] - TT[-2])/(tt[-1] - tt[-2])
        delay=tt[np.argmax(dTdt)]

        return delay
#        gas=ct.Solution('mechanism')
#        gas.TPX = T,P,X
#        r = ct.IdealGasConstPressureReactor(gas)
#        sim = ct.ReactorNet([r])
#    
#        tt = []
#        TT = []
#        t = 0.0
#        
#        Rmax = np.zeros(gas.n_reactions)
#        while t < 0.05:
#            t = sim.step()
#            tt.append(1000 * t)
#            TT.append(r.T)
#            rnet = abs(gas.net_rates_of_progress)
#            rnet /= max(rnet)
#            Rmax = np.maximum(Rmax, rnet)
#            
#        delay=[]
#        dTdt=np.zeros(np.array(tt).shape,np.float)
#        dTdt[0:-1]=np.diff(TT)/np.diff(tt)
#        dTdt[-1]=(TT[-1] - TT[-2])/(tt[-1] - tt[-2])
#        
#        delay=tt[np.argmax(dTdt)]
#        return delay
    elif options=='cv':
        #gas.TPX = T,P,X
        #r = ct.IdealGasReactor(gas)
        #sim = ct.ReactorNet([r])
        tt = []
        TT = []
        results=stcv.ShockTube(mechanism,['N2'],P,T,X,0,10000.0,'adiabatic',reactorType='cv')
        #t = 0.0
        #while t < 0.02:
            #t = sim.step()
            #tt.append(1000 * t)
            #TT.append(r.T)
        delay=[]
        tt=results.solution['time'].values
        TT=results.solution['temperature'].values
        PP=results.solution['pressure'].values
        dTdt=np.zeros(np.array(tt).shape,np.float)
        dTdt[0:-1]=np.diff(TT)/np.diff(tt)
        dTdt[-1]=(TT[-1] - TT[-2])/(tt[-1] - tt[-2])
        delay=tt[np.argmax(dTdt)]
        dPdt=np.zeros(np.array(tt).shape,np.float)
        dPdt[0:-1]=np.diff(PP)/np.diff(tt)
        dPdt[-1]=(PP[-1] - PP[-2])/(tt[-1] - tt[-2])
        delayp=tt[np.argmax(dPdt)]
        if igdelay_type=='T':
            return delay
        elif igdelay_type=='P':
            return delayp
        
        
        
        
        