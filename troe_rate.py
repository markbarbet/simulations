# -*- coding: utf-8 -*-
"""
Created on Fri Feb 02 18:12:43 2018

@author: Mark Barbet

"""


import cantera as ct
import numpy as np

def rate(kl,kh,a,t3,t1,t2,Temperature,pressure,params=[]):
    Fcent=(1-a)*np.exp(-Temperature/t3)+a*np.exp(-Temperature/t1)+np.exp(-t2/Temperature)
    
    c = -0.4-0.67*np.log10(Fcent)
    n = 0.75 - 1.27 *np.log10(Fcent)
    d = 0.14
    
    Pr = kl*pressure/ct.gas_constant/Temperature/kh
    
    F = 10.0**((((1.0+(np.divide(np.log10(Pr)+c,n-d*(np.log10(Pr)+c)))**2.0))**(-1.0))*np.log10(Fcent))
    
    k = kh*Pr*F/(1+Pr)
    if 'Pr' in params:
        print(Pr)
    return k

    