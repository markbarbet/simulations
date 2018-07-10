# -*- coding: utf-8 -*-
"""
Created on Sat Feb 03 16:15:41 2018

@author: Mark Barbet
"""
import numpy as np


def rate(A,n,Ea,T):
    #calculates the rate constant using calories in activation energy
    
    rate=A*(T**n)*np.exp(-Ea/T/1.9858775)
    return rate